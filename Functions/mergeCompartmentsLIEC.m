function mergedModel = mergeCompartments(model, mergeRules)
% Merges specified compartments in a COBRA model
%
% INPUT:
%   model       - COBRA model structure
%   mergeRules  - nx2 cell array {oldComp, newComp}
%                 newComp will be the output compartment
%
% OUTPUT:
%   mergedModel - model with merged compartments
%
% EXAMPLE:
%   rules = {'[lu]',  '[luLI]';
%            '[luC]', '[luLI]';
%            '[bpC]', '[e]'};
%   mergedModel = mergeCompartments(ColonModelFemale, rules)
%
% by Sena Erimli, 2026

mergedModel = model;
nOrigMets = length(model.mets);
nOrigRxns = length(model.rxns);

%% STEP 1: Validate input
if size(mergeRules, 2) ~= 2
    error('mergeRules must be an nx2 cell array: {oldComp, newComp}')
end
fprintf('=== Merge rules ===\n')
for r = 1:size(mergeRules, 1)
    fprintf('  %s -> %s\n', mergeRules{r,1}, mergeRules{r,2})
end

%% STEP 2: Rename metabolites
fprintf('\n=== Renaming metabolites ===\n')
for r = 1:size(mergeRules, 1)
    oldComp = mergeRules{r, 1};
    newComp = mergeRules{r, 2};

    oldIdx = find(contains(mergedModel.mets, oldComp));
    if isempty(oldIdx)
        fprintf('  No metabolites found with %s\n', oldComp)
        continue
    end

    mergedModel.mets(oldIdx) = strrep(mergedModel.mets(oldIdx), oldComp, newComp);

    allFields = fieldnames(mergedModel);
    for i = 1:length(allFields)
        field = allFields{i};
        val   = mergedModel.(field);
        if iscell(val) && size(val,1) == nOrigMets && ...
                all(cellfun(@ischar, val))
            mergedModel.(field)(oldIdx) = strrep(...
                mergedModel.(field)(oldIdx), oldComp, newComp);
        end
    end

    fprintf('  %s -> %s : %d metabolites renamed\n', ...
        oldComp, newComp, length(oldIdx))
end

%% STEP 3: Rename reactions
fprintf('\n=== Renaming reactions ===\n')
for r = 1:size(mergeRules, 1)
    oldComp = mergeRules{r, 1};
    newComp = mergeRules{r, 2};

    oldRxnIdx = find(contains(mergedModel.rxns, oldComp));
    if isempty(oldRxnIdx)
        continue
    end

    allFields = fieldnames(mergedModel);
    for i = 1:length(allFields)
        field = allFields{i};
        val   = mergedModel.(field);
        if iscell(val) && size(val,1) == nOrigRxns && ...
                all(cellfun(@ischar, val))
            mergedModel.(field)(oldRxnIdx) = strrep(...
                mergedModel.(field)(oldRxnIdx), oldComp, newComp);
        end
    end

    fprintf('  %s -> %s : %d reactions renamed\n', ...
        oldComp, newComp, length(oldRxnIdx))
end



%% STEP 4: Merge duplicate metabolites
fprintf('\n=== Merging duplicate metabolites ===\n')
[~, ia_met, ic_met] = unique(mergedModel.mets, 'stable');
dupMetCount = nOrigMets - length(ia_met);
fprintf('Duplicate metabolites to merge: %d\n', dupMetCount)

if dupMetCount > 0
    newS = zeros(length(ia_met), size(mergedModel.S, 2));
    for i = 1:length(ia_met)
        dupIdx = find(ic_met == i);
        newS(i, :) = sum(mergedModel.S(dupIdx, :), 1);
    end
    mergedModel.S = sparse(newS);

    allFields = fieldnames(mergedModel);
    for i = 1:length(allFields)
        field = allFields{i};
        if strcmp(field, 'S')  || strcmp(field, 'C'), continue, end
        val = mergedModel.(field);
        if (iscell(val) || isnumeric(val)) && size(val,1) == nOrigMets
            mergedModel.(field) = val(ia_met, :);
        end
    end

    fprintf('Metabolites before: %d | after: %d\n', ...
        nOrigMets, length(mergedModel.mets))
end

%% STEP 4: Remove self-loop reactions
fprintf('\n=== Removing self-loop reactions ===\n')
selfLoopIdx = false(length(mergedModel.rxns), 1);

for i = 1:length(mergedModel.rxns)
    rxnID = mergedModel.rxns{i};

    % Skip demand, sink, and exchange reactions
    % These legitimately have metabolites in one compartment
    if startsWith(rxnID, 'DM_')   || ...
       startsWith(rxnID, 'sink_') || ...
       startsWith(rxnID, 'EX_')
        continue
    end

    % Use S matrix directly to find participating metabolites
    metIdx = find(mergedModel.S(:, i) ~= 0);
    if isempty(metIdx)
        continue
    end

    % Extract compartments from metabolite names
    rxnMets = mergedModel.mets(metIdx);
    comps   = regexp(rxnMets, '\[(\w+)\]', 'tokens');
    comps   = cellfun(@(x) x{1}{1}, comps, 'UniformOutput', false);

    % Self-loop: transport reaction where all metabolites
    % are now in the same compartment after renaming
    if length(unique(comps)) == 1 && length(metIdx) > 1
        selfLoopIdx(i) = true;
        fprintf('  Removing self-loop: %s\n', rxnID)
    end
end

fprintf('Self-loop reactions removed: %d\n', sum(selfLoopIdx))

%% STEP 6: Remove duplicate reactions
fprintf('\n=== Removing duplicate reactions ===\n')
[~, ia_rxn, ic_rxn] = unique(mergedModel.rxns, 'stable');
dupRxnCount = nOrigRxns - length(ia_rxn);
fprintf('Duplicate reactions to remove: %d\n', dupRxnCount)

if dupRxnCount > 0
    for i = 1:length(ia_rxn)
        dupIdx = find(ic_rxn == i);
        if length(dupIdx) > 1
            %Stoichiometry check
            ref = full(mergedModel.S(:, dupIdx(1)));
            for j = 2:length(dupIdx)
                if ~isequal(ref, full(mergedModel.S(:, dupIdx(j))))
                    fprintf('  WARNING: Different stoich for duplicate: %s\n', ...
                        mergedModel.rxns{ia_rxn(i)})
                end
            end
            % --- Gene rules ---
            if isfield(mergedModel, 'grRules')
                rules = mergedModel.grRules(dupIdx);
                rules = rules(~cellfun(@isempty, rules));
                rules = unique(rules);

                if isempty(rules)
                    mergedModel.grRules{dupIdx(1)} = '';
                elseif length(rules) == 1
                    mergedModel.grRules{dupIdx(1)} = rules{1};
                else
                    wrapped = cellfun(@(r) ['(' r ')'], rules, ...
                                      'UniformOutput', false);
                    mergedModel.grRules{dupIdx(1)} = strjoin(wrapped, ' or ');
                    fprintf('  Gene rules merged for: %s\n', ...
                        mergedModel.rxns{ia_rxn(i)})
                end
            end         
        end
    end

    % Keep first occurrence - gene rules already updated above
    mergedModel.S = mergedModel.S(:, ia_rxn);
    
    if isfield(mergedModel, 'C') && size(mergedModel.C, 2) == nOrigRxns
        mergedModel.C = mergedModel.C(:, ia_rxn);
    end

    allFields = fieldnames(mergedModel);
    for i = 1:length(allFields)
        field = allFields{i};
        if strcmp(field, 'S') || strcmp(field, 'C'), continue, end
        val = mergedModel.(field);
        if (iscell(val) || isnumeric(val)) && size(val,1) == nOrigRxns
            mergedModel.(field) = val(ia_rxn, :);
        elseif isnumeric(val) && size(val,2) == nOrigRxns
            mergedModel.(field) = val(:, ia_rxn);
        end
    end

    fprintf('Reactions before: %d | after: %d\n', ...
        nOrigRxns, length(mergedModel.rxns))
end

%% STEP 7: Update comps and compNames fields
compFields = {'comps', 'compNames'};
for f = 1:length(compFields)
    field = compFields{f};
    if isfield(mergedModel, field)
        for r = 1:size(mergeRules, 1)
            oldComp = strrep(mergeRules{r,1}, '[', '');
            oldComp = strrep(oldComp, ']', '');
            newComp = strrep(mergeRules{r,2}, '[', '');
            newComp = strrep(newComp, ']', '');
            mergedModel.(field) = strrep(mergedModel.(field), oldComp, newComp);
        end
        mergedModel.(field) = unique(mergedModel.(field), 'stable');
    end
end

%% STEP 8: Verify
sol = optimizeCbModel(mergedModel);
fprintf('\n=== Verification ===\n')
fprintf('Biomass after merge: %.4f\n', sol.f)
fprintf('Remaining compartments:\n')
met_comps = regexp(mergedModel.mets, '\[(\w+)\]', 'tokens');
met_comps = met_comps(~cellfun(@isempty, met_comps));
met_comps = cellfun(@(x) x{1}{1}, met_comps, 'UniformOutput', false);
disp(unique(met_comps))
end