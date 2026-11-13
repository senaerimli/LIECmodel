function [boundDiffs, idMap, uniqueRxns, uniqueMets, updatedModel1] = compareStoichiometry(model1, model2, name1, name2, updateFlag)
% COMPARESTOICHIOMETRY Compares models based on the S-matrix, identifies boundary
% differences, conditionally updates model1, and extracts categorized metabolites.
%
% INPUTS:
%   model1, model2 : COBRA model structures
%   name1, name2   : (Optional) Strings for table headers
%   updateFlag     : (Optional) Set to 1 to automatically update model1's bounds.
%
% OUTPUTS:
%   boundDiffs     : Table of mathematically identical reactions with DIFFERENT bounds.
%   idMap          : Master dictionary translating Model 1 to Model 2.
%   uniqueRxns     : Struct of tables for reactions unique to only one model.
%   updatedModel1  : Model 1 (updated if updateFlag=1, otherwise original).
%   uniqueMets     : Nested struct containing lists of unique base metabolites 
%                    organized specifically by their Category pathways!

    if nargin < 5 || isempty(updateFlag)
        updateFlag = 0;
    end
    if nargin < 4 || isempty(name2)
        name2 = 'Model2';
    end
    if nargin < 3 || isempty(name1)
        name1 = 'Model1';
    end

    fprintf('\n======================================================\n');
    fprintf(' Running Flag-Controlled S-Matrix Comparison...\n');
    fprintf(' %s vs %s\n', name1, name2);
    fprintf('======================================================\n');

    updatedModel1 = model1;

    % 1. Normalize Extreme Bounds
    updatedModel1.lb(updatedModel1.lb == -1000000) = -1000;
    updatedModel1.ub(updatedModel1.ub ==  1000000) =  1000;
    model2.lb(model2.lb == -1000000) = -1000;
    model2.ub(model2.ub ==  1000000) =  1000;

    % 2. Generate Stoichiometric Signatures
    fprintf('Extracting signatures and compartments for %s...\n', name1);
    sigs1 = getRxnSignatures(updatedModel1);
    
    fprintf('Extracting signatures and compartments for %s...\n', name2);
    sigs2 = getRxnSignatures(model2);

    % 3. Find Exact Mathematical Intersections
    [commonSigs, idx1, idx2] = intersect(sigs1, sigs2);
    fprintf('\nFound %d mathematically identical reactions.\n', length(commonSigs));

    % 4. Categorize the shared reactions
    sharedCategories = getRxnCategories(updatedModel1, idx1);

    % 5. Build the Master ID Translation Map (idMap)
    ID_in_Model1 = updatedModel1.rxns(idx1);
    ID_in_Model2 = model2.rxns(idx2);
    Formula = printRxnFormula(updatedModel1, 'rxnAbbrList', ID_in_Model1, 'printFlag', false);
    
    LB_1 = updatedModel1.lb(idx1);
    UB_1 = updatedModel1.ub(idx1);
    LB_2 = model2.lb(idx2);
    UB_2 = model2.ub(idx2);
    
    idMap = table(sharedCategories, ID_in_Model1, ID_in_Model2, Formula, LB_1, UB_1, LB_2, UB_2, ...
        'VariableNames', {'Category', ['ID_', name1], ['ID_', name2], 'Formula', ...
                          ['LB_', name1], ['UB_', name1], ...
                          ['LB_', name2], ['UB_', name2]});
                          
    idMap = sortrows(idMap, 'Category');

    % 6. Filter for Boundary Differences
    diff_mask = (idMap.(['LB_', name1]) ~= idMap.(['LB_', name2])) | ...
                (idMap.(['UB_', name1]) ~= idMap.(['UB_', name2]));
    boundDiffs = idMap(diff_mask, :);

    if isempty(boundDiffs)
        fprintf('Perfect Match! All shared reactions already have identical boundaries.\n');
    else
        fprintf('Warning: Found %d reactions with DIFFERENT bounds.\n', height(boundDiffs));
        
        if updateFlag == 1
            fprintf(' -> updateFlag=1 detected. Syncing %s bounds to match %s...\n', name1, name2);
            for i = 1:height(boundDiffs)
                rxn_ID_to_update = boundDiffs.(['ID_', name1]){i};
                new_lb = boundDiffs.(['LB_', name2])(i);
                new_ub = boundDiffs.(['UB_', name2])(i);
                
                rxn_loc = strcmp(updatedModel1.rxns, rxn_ID_to_update);
                updatedModel1.lb(rxn_loc) = new_lb;
                updatedModel1.ub(rxn_loc) = new_ub;
            end
            fprintf(' -> Successfully generated patched model.\n');
        else
            fprintf(' -> updateFlag=0 (or not set). Model 1 remains unchanged.\n');
        end
    end

    % 7. Categorize Truly Unique Reactions
    unmatched1 = setdiff(1:length(updatedModel1.rxns), idx1);
    unmatched2 = setdiff(1:length(model2.rxns), idx2);
    
    uniqueRxns = struct();
    if ~isempty(unmatched1)
        cat1 = getRxnCategories(updatedModel1, unmatched1);
        form1 = printRxnFormula(updatedModel1, 'rxnAbbrList', updatedModel1.rxns(unmatched1), 'printFlag', false);
        uniqueRxns.(name1) = sortrows(table(cat1, updatedModel1.rxns(unmatched1), form1, updatedModel1.lb(unmatched1), updatedModel1.ub(unmatched1), ...
            'VariableNames', {'Category', 'Reaction_ID', 'Formula', 'lb', 'ub'}), 'Category');
    else, uniqueRxns.(name1) = table(); end

    if ~isempty(unmatched2)
        cat2 = getRxnCategories(model2, unmatched2);
        form2 = printRxnFormula(model2, 'rxnAbbrList', model2.rxns(unmatched2), 'printFlag', false);
        uniqueRxns.(name2) = sortrows(table(cat2, model2.rxns(unmatched2), form2, model2.lb(unmatched2), model2.ub(unmatched2), ...
            'VariableNames', {'Category', 'Reaction_ID', 'Formula', 'lb', 'ub'}), 'Category');
    else, uniqueRxns.(name2) = table(); end

    % --- NEW: 8. Extract Categorized Base Metabolites ---
    uniqueMets = struct();
    
    % Get metabolites for boundary differences, grouped by category
    uniqueMets.BoundDifferences = getCategorizedMets(updatedModel1, boundDiffs, ['ID_', name1]);
    
    % Get metabolites for unique reactions, grouped by category
    uniqueMets.(name1) = getCategorizedMets(updatedModel1, uniqueRxns.(name1), 'Reaction_ID');
    uniqueMets.(name2) = getCategorizedMets(model2, uniqueRxns.(name2), 'Reaction_ID');

    fprintf('======================================================\n');
end

% =========================================================================
% HELPER 1: Generates a unique string for a reaction's stoichiometry
% =========================================================================
function sigs = getRxnSignatures(model)
    sigs = cell(length(model.rxns), 1);
    for i = 1:length(model.rxns)
        metsIdx = find(model.S(:, i));
        if isempty(metsIdx), sigs{i} = 'EMPTY_RXN'; continue; end
        mets = model.mets(metsIdx);
        coeffs = full(model.S(metsIdx, i));
        [sortedMets, sortIdx] = sort(mets);
        sortedCoeffs = coeffs(sortIdx);
        
        sig = '';
        for j = 1:length(sortedMets)
            sig = [sig, sprintf('%s:%.5g|', sortedMets{j}, sortedCoeffs(j))];
        end
        sigs{i} = sig;
    end
end

% =========================================================================
% HELPER 2: Categorizes reactions into Exchange, Transport, or Internal
% =========================================================================
function categories = getRxnCategories(model, rxnIndices)
    metComps = cell(length(model.mets), 1);
    for i = 1:length(model.mets)
        tok = regexp(model.mets{i}, '\[([a-zA-Z0-9_]+)\]$', 'tokens');
        if isempty(tok), tok = regexp(model.mets{i}, '_([a-zA-Z0-9]+)$', 'tokens'); end
        if ~isempty(tok), metComps{i} = tok{1}{1}; else metComps{i} = 'unknown'; end
    end

    categories = cell(length(rxnIndices), 1);
    for k = 1:length(rxnIndices)
        metIdx = find(model.S(:, rxnIndices(k)));
        if isempty(metIdx), categories{k} = 'Empty_Reaction'; continue; end
        
        rxnComps = unique(metComps(metIdx));
        if length(metIdx) == 1, categories{k} = sprintf('Exchange_from_%s', rxnComps{1});
        elseif length(rxnComps) == 1, categories{k} = sprintf('Internal_%s', rxnComps{1});
        elseif length(rxnComps) == 2, sortedComps = sort(rxnComps); categories{k} = sprintf('%s_AND_%s', sortedComps{1}, sortedComps{2});
        else, categories{k} = 'Multi_Compartment'; end
    end
end

% =========================================================================
% HELPER 3: Extracts base metabolites and groups them by Category
% =========================================================================
function catMets = getCategorizedMets(model, rxnTable, idColName)
    catMets = struct();
    if isempty(rxnTable)
        return;
    end
    
    % Find all unique categories in the provided table
    uniqueCats = unique(rxnTable.Category);
    
    for i = 1:length(uniqueCats)
        currentCat = uniqueCats{i};
        
        % Ensure the category name is a valid MATLAB struct field name
        cleanCatName = regexprep(currentCat, '[^a-zA-Z0-9_]', '');
        
        % Isolate the reactions belonging to this specific category
        subTable = rxnTable(strcmp(rxnTable.Category, currentCat), :);
        rxnIDs = subTable.(idColName);
        
        baseMets = {};
        for j = 1:length(rxnIDs)
            rxnIdx = find(strcmp(model.rxns, rxnIDs{j}));
            if ~isempty(rxnIdx)
                metIdx = find(model.S(:, rxnIdx));
                rxnMets = model.mets(metIdx);
                
                % Strip compartment tags
                cleanedMets = regexprep(rxnMets, '\[([a-zA-Z0-9_]+)\]$|_([a-zA-Z0-9]+)$', '');
                baseMets = [baseMets; cleanedMets];
            end
        end
        % Save the unique list of metabolites directly into the category field
        catMets.(cleanCatName) = unique(baseMets);
    end
end