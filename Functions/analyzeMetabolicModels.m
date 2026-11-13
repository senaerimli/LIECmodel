function [Results, UpdatedModels] = analyzeMetabolicModels(modelsStruct, objectiveName, resultsPath, setName)
% analyzeMetabolicModels Runs FBA, shadow prices, dead ends, and sanity checks on multiple models.
%
% INPUTS:
%   modelsStruct  - A struct where each field is a COBRA model (e.g., modelsStruct.ColonModelFemale = model;)
%   objectiveName - String specifying the objective to set (e.g., 'biomass_maintenance')
%   resultsPath   - String specifying the directory to save the Excel files.
%   setName       - Directory of the output excel files
% OUTPUTS:
%   Results       - Struct containing all outputs (deadEnds, FBA solutions, shadowPrices, sanityChecks)
%   UpdatedModels - Struct containing the models after the objective has been changed.

% 1. Initialization
    modelNames = fieldnames(modelsStruct);
    numModels = length(modelNames);
    
    Results = struct();
    UpdatedModels = struct();
    
    fprintf('\n=== Starting Metabolic Model Analysis ===\n');
    fprintf('Processing %d models...\n', numModels);
    
    % 2. Process each model (FBA, Dead Ends, Sanity, Shadow Prices)
    for i = 1:numModels
        currentName = modelNames{i};
        currentModel = modelsStruct.(currentName);
        
        fprintf('\n--- Analyzing [%d/%d]: %s ---\n', i, numModels, currentName);
        
        % Detect Dead Ends
        deadIdx = detectDeadEnds(currentModel);
        Results.deadEnds.(currentName) = deadIdx;
        fprintf('Total Dead Ends: %d\n', length(deadIdx));
        
        % Change Objective
        currentModel = changeObjective(currentModel, objectiveName);
        UpdatedModels.(currentName) = currentModel;
        
        % Run FBA
        try
            sol = optimizeCbModel(currentModel);
            Results.fba.(currentName) = sol;
            fprintf('FBA Objective Value: %.4f\n', sol.f);
        catch ME
            fprintf('Warning: FBA failed for %s. Error: %s\n', currentName, ME.message);
            Results.fba.(currentName) = [];
        end
        
        % Calculate Shadow Prices
        try
            [more, zero, less] = calculateShadowPrice(currentModel);
            Results.shadowPrices.moreThanOne.(currentName) = more;
            Results.shadowPrices.equalsToZero.(currentName) = zero;
            Results.shadowPrices.lessThanOne.(currentName) = less;
        catch ME
            fprintf('Warning: Shadow price calculation failed for %s.\n', currentName);
        end
        
        % Perform Sanity Checks
       try
            % Generate the safe sheet name (e.g., 'Model_001')
            safeSheetName = sprintf('Model_%03d', i);
            
            % Pass it as the 3rd argument to your updated function
            sanity = sanityChecksFromTutorial(currentModel, resultsPath, safeSheetName);
            Results.sanityChecks.(currentName) = sanity;
        catch ME
            fprintf('⚠️ Warning: Sanity checks failed for %s.\n   -> ERROR: %s\n', currentName, ME.message);
        end
    end
    
    % 3. Export Results to Excel with Safe "Index" Mapping
    fprintf('\n=== Starting Excel Export ===\n');
    
% Use sprintf to attach the setName to the beginning of the Excel files
    sanityExcel     = fullfile(resultsPath, sprintf('%s_SanityChecks_Results.xlsx', setName));
    shadowMoreExcel = fullfile(resultsPath, sprintf('%s_ShadowPrices_MoreThanOne.xlsx', setName));
    shadowZeroExcel = fullfile(resultsPath, sprintf('%s_ShadowPrices_EqualsToZero.xlsx', setName));
    shadowLessExcel = fullfile(resultsPath, sprintf('%s_ShadowPrices_LessThanOne.xlsx', setName));
    
    % Create a mapping array for the Master Index Sheet
    sheetMapping = cell(numModels, 2);
    
    for i = 1:numModels
        currentName = modelNames{i};
        
        % Create a perfectly safe, short sheet name (e.g., 'Model_001')
        safeSheetName = sprintf('Model_%03d', i); 
        
        % Save to our mapping array
        sheetMapping{i, 1} = safeSheetName;
        sheetMapping{i, 2} = currentName;
        
        fprintf('Exporting %s as %s...\n', currentName, safeSheetName);
        
        % Export Sanity Checks
        if isfield(Results, 'sanityChecks') && isfield(Results.sanityChecks, currentName) && ~isempty(Results.sanityChecks.(currentName))
            T_sanity = cell2table(Results.sanityChecks.(currentName), 'VariableNames', {'Test', 'Result'});
            writetable(T_sanity, sanityExcel, 'Sheet', safeSheetName);
        end
        
        % Export Shadow Prices using helper function
        exportShadowPrices(Results.shadowPrices.moreThanOne, currentName, shadowMoreExcel, safeSheetName);
        exportShadowPrices(Results.shadowPrices.equalsToZero, currentName, shadowZeroExcel, safeSheetName);
        exportShadowPrices(Results.shadowPrices.lessThanOne, currentName, shadowLessExcel, safeSheetName);
    end
    
    % 4. Write the Master Index Sheet to all files
    fprintf('\nWriting 00_INDEX master mapping sheets...\n');
    T_Index = cell2table(sheetMapping, 'VariableNames', {'Excel_Sheet_Name', 'Full_MATLAB_Model_Name'});
    
    % Suppress warnings about adding new sheets to Excel
    warning('off', 'MATLAB:xlswrite:AddSheet');
    
    writetable(T_Index, sanityExcel, 'Sheet', '00_INDEX');
    writetable(T_Index, shadowMoreExcel, 'Sheet', '00_INDEX');
    writetable(T_Index, shadowZeroExcel, 'Sheet', '00_INDEX');
    writetable(T_Index, shadowLessExcel, 'Sheet', '00_INDEX');
    
    warning('on', 'MATLAB:xlswrite:AddSheet');
    
    fprintf('✅ Excel export completed successfully!\n');
end
