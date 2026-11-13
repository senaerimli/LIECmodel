function model = updateModelBounds(model, updateTable, targetColName)
% UPDATEMODELBOUNDS Updates model boundaries based on a reference table.
% 
% model: The COBRA model struct to be updated.
% updateTable: The table containing new boundary values.
% targetColName: The string name of the column in updateTable to use for the update.

    % 1. Find the overlap between the table and the model
    [isMatch, modelIdx] = ismember(updateTable.Reaction, model.rxns);
    
    % 2. Safety check: ensure we found matches
    matchCount = sum(isMatch);
    if matchCount == 0
        warning('No matching reactions found between table and model.');
        return;
    end
    
    % 3. Extract the new values from the table
    % Only take rows that actually exist in the model
    newValues = updateTable.(targetColName)(isMatch);
    
    % 4. Apply the update
    % Use the corresponding indices in the model
    validModelIdx = modelIdx(isMatch);
    
    % Check if we are updating UB or LB based on column name hints
    if contains(targetColName, 'ub', 'IgnoreCase', true)
        model.ub(validModelIdx) = newValues;
        fprintf('Updated Upper Bounds for %d reactions.\n', matchCount);
    elseif contains(targetColName, 'lb', 'IgnoreCase', true)
        model.lb(validModelIdx) = newValues;
        fprintf('Updated Lower Bounds for %d reactions.\n', matchCount);
    else
        % Default: If the column name isn't clear, you might need to specify 
        % but for your case, we assume it matches the column you passed.
        model.ub(validModelIdx) = newValues; 
    end
end