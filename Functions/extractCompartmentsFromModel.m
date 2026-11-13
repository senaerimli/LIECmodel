function subModel = extractCompartmentsFromModel(model, compartments, prefixes, stripPrefix)
% Extracts reactions involving specified compartments and/or name prefixes
%
% INPUT:
%   model        - source COBRA model (e.g., WBM)
%   compartments - cell array of compartment tags, or {} to skip
%                  e.g., {'[luLI]', '[fe]', '[bpC]', '[d]', '[luC]'}
%   prefixes     - (optional) cell array of prefixes to match in rxn/met names
%                  e.g., {'Colon_'} or {} to skip
%                  if not provided, defaults to {'Colon_'}
%
%   optional-
%   stripPrefix  - if "true" the prefix will be striped from rxns and mets
%
% OUTPUT:
%   subModel     - extracted model structure with all fields preserved
%
% EXAMPLES:
%   % Original behaviour - compartments + Colon_ prefix
%   subModel = extractCompartmentsFromModel(modelWBM, ...
%       {'[luLI]', '[fe]', '[bpC]', '[d]', '[luC]'})
%
%   % Compartments only - no prefix matching
%   subModel = extractCompartmentsFromModel(modelWBM, ...
%       {'[luLI]', '[fe]'}, {})
%
%   % Custom prefixes
%   subModel = extractCompartmentsFromModel(modelWBM, ...
%       {'[luLI]'}, {'Colon_', 'LI_'})
%
%   % Prefixes only - no compartment matching
%   subModel = extractCompartmentsFromModel(modelWBM, ...
%       {}, {'Colon_'})
%
%   compartments = {'[luLI]', '[fe]', '[bpC]', '[d]', '[luC]'};
%   prefixes = {'Colon_'};
%   subModel = extractCompartmentsFromModel(model, compartments, prefixes, true);
%
% by Sena Erimli, 2026

%% STEP 1: Validate and set defaults
if nargin < 2
    error('At least model and compartments must be provided.')
end
if nargin < 3
    prefixes = {'Colon_'};  % default behaviour
end

% Both empty is an error
if isempty(compartments) && isempty(prefixes)
    error('At least one of compartments or prefixes must be non-empty.')
end

% Print what we are doing
fprintf('=== extractCompartmentsFromModel ===\n')
if ~isempty(compartments)
    fprintf('Compartments: %s\n', strjoin(compartments, ', '))
else
    fprintf('Compartments: none\n')
end
if ~isempty(prefixes)
    fprintf('Prefixes:     %s\n', strjoin(prefixes, ', '))
else
    fprintf('Prefixes:     none\n')
end

%% STEP 2: Build match index
matchIdx = false(length(model.rxns), 1);

% Match reaction names by prefix
for p = 1:length(prefixes)
    matchIdx = matchIdx | startsWith(model.rxns, prefixes{p});
end

% Match reactions involving target compartments
for c = 1:length(compartments)
    metIdx = find(contains(model.mets, compartments{c}));
    if ~isempty(metIdx)
        rxnIdx = any(model.S(metIdx, :) ~= 0, 1)';
        matchIdx = matchIdx | rxnIdx;
    end
end

% Match reactions involving metabolites with target prefixes
for p = 1:length(prefixes)
    prefMetIdx = find(startsWith(model.mets, prefixes{p}));
    if ~isempty(prefMetIdx)
        rxnIdx = any(model.S(prefMetIdx, :) ~= 0, 1)';
        matchIdx = matchIdx | rxnIdx;
    end
end

%% STEP 3: Extract indices
rxnIdx = find(matchIdx);
metIdx = find(any(model.S(:, rxnIdx) ~= 0, 2));

if isempty(rxnIdx)
    error('No reactions matched. Check compartment tags and prefixes.')
end

%% STEP 4: Build subModel
subModel = struct();

% Core fields handled explicitly
subModel.S  = model.S(metIdx, rxnIdx);
subModel.lb = model.lb(rxnIdx);
subModel.ub = model.ub(rxnIdx);
subModel.c  = model.c(rxnIdx);
subModel.b  = model.b(metIdx);

skipFields = {'S', 'lb', 'ub', 'c', 'b'};

allFields = fieldnames(model);
for i = 1:length(allFields)
    field = allFields{i};

    if ismember(field, skipFields), continue, end

    val = model.(field);

    % Non-array fields - copy as-is
    if ~iscell(val) && ~isnumeric(val)
        subModel.(field) = val;
        continue
    end

    nRows = size(val, 1);
    nCols = size(val, 2);

    if nRows == length(model.rxns) && nCols == length(model.mets)
        subModel.(field) = val(rxnIdx, metIdx);
    elseif nRows == length(model.mets) && nCols == length(model.rxns)
        subModel.(field) = val(metIdx, rxnIdx);
    elseif nRows == length(model.rxns)
        subModel.(field) = val(rxnIdx, :);
    elseif nRows == length(model.mets)
        subModel.(field) = val(metIdx, :);
    elseif nCols == length(model.rxns)
        subModel.(field) = val(:, rxnIdx);
    elseif nCols == length(model.mets)
        subModel.(field) = val(:, metIdx);
    else
        subModel.(field) = val;
    end
end

% Strip prefixes from extracted model if requested
if nargin >= 4 && stripPrefix
    for p = 1:length(prefixes)
        pattern = ['^' prefixes{p}];
        subModel.rxns = regexprep(subModel.rxns, pattern, '');
        subModel.mets = regexprep(subModel.mets, pattern, '');
    end
    fprintf('Prefix(es) stripped from reaction and metabolite names\n')
end

fprintf('Extracted reactions:   %d\n', length(subModel.rxns))
fprintf('Extracted metabolites: %d\n', length(subModel.mets))
fprintf('Fields preserved:      %d\n', length(fieldnames(subModel)))
end