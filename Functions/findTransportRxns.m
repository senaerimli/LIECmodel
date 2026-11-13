function transportRxnsTable = findTransportRxns(model, comp1, comp2, filterCustomBounds)
% findTransportRxns Identifies reactions that transport metabolites 
% between two specified compartments.
%
% INPUTS:
%   model              - COBRA model structure
%   comp1              - String tag for the first compartment (e.g., '[luLI]')
%   comp2              - String tag for the second compartment (e.g., '[fe]')
%   filterCustomBounds - (Optional) Boolean. If true, only returns reactions 
%                        where lb or ub are NOT -1000, 0, or 1000. 
%                        Default is true.
%
% OUTPUT:
%   transportRxnsTable - A MATLAB table with the reaction ID, formula, and bounds.
% USAGE
%   transRxnsT = findTransportRxns(model, '[luLI]', '[fe]', true)


% Default to filtering for custom bounds if not specified
if nargin < 4
    filterCustomBounds = true;
end

model.lb(model.lb == -1000000) = -1000;
model.ub(model.ub ==  1000000) =  1000;

standard_bounds = [-1000, 0, 1000];
results_cell = {};

for i = 1:length(model.rxns)
    % 1. Identify metabolites involved in this specific reaction
    met_indices = find(model.S(:, i));
    rxn_mets = model.mets(met_indices);

    % 2. Check if the reaction bridges both compartments
    has_comp1 = any(contains(rxn_mets, comp1));
    has_comp2 = any(contains(rxn_mets, comp2));

    if has_comp1 && has_comp2
        lb_val = model.lb(i);
        ub_val = model.ub(i);

        % 3. Apply bounds filter if requested
        if filterCustomBounds
            is_lb_standard = ismember(lb_val, standard_bounds);
            is_ub_standard = ismember(ub_val, standard_bounds);

            % If both bounds are standard (-1000, 0, 1000), skip to the next reaction
            if is_lb_standard && is_ub_standard
                continue;
            end
        end

        % 4. Extract data for matching reactions
        rxn_id = model.rxns{i};
        
        % Generate the readable reaction formula
        formula_cell = printRxnFormula(model, 'rxnAbbrList', {rxn_id}, 'printFlag', false);
        formula = formula_cell{1};

        % Store the results
        results_cell(end+1, :) = {rxn_id, formula, lb_val, ub_val};
    end
end

% 5. Format and return the output as a table
if isempty(results_cell)
    fprintf('No matching transport reactions found between %s and %s.\n', comp1, comp2);
    transportRxnsTable = table(); % Return an empty table
else
    transportRxnsTable = cell2table(results_cell, ...
        'VariableNames', {'Reaction_ID', 'Formula', 'lb', 'ub'});
    fprintf('Found %d transport reactions.\n', height(transportRxnsTable));
end
end