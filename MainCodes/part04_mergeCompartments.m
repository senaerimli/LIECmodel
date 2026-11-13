%%Merge Compartments
%by Sena Erimli

clear all
close all 
clc

initCobraToolbox(false)
addpath(genpath('.'));

% 1. DYNAMIC PATHS
scriptPath = mfilename('fullpath');

% If running via F9 or Live Script, mfilename is empty. Grab the active file from the Editor instead.
if isempty(scriptPath)
    scriptPath = matlab.desktop.editor.getActiveFilename;
end

% Calling fileparts twice steps up
basePath = fileparts(fileparts(scriptPath));

global modelsPath resultsPath 
modelsPath  = fullfile(basePath, 'Models');
resultsPath = fullfile(basePath, 'Results');

addpath(modelsPath, resultsPath);

if 1
    %LPSolver = 'tomlab_cplex';
    %QPSolver = 'tomlab_cplex';
    LPSolver = 'gurobi';
    QPSolver = 'gurobi';
else
    LPSolver = 'ibm_cplex';
    QPSolver = 'ibm_cplex';
end

[solverOK, solverInstalled] = changeCobraSolver(LPSolver, 'LP',0,1);
[solverOK, solverInstalled] = changeCobraSolver(QPSolver, 'QP',0,1);

%% Load the necessary models
load ColonModelMale.mat;
load ColonModelFemale.mat;

load ColonModelMaleupd.mat;
load ColonModelFemaleupd.mat;

load male_ColonModelDietFe_noConst   %colonModelMale_withDietFe
load female_ColonModelDietFe_noConst %colonModelFemale_withDietFe

load colonModelMale_withDietFeupd
load colonModelFemale_withDietFeupd

%% Merge boundary compartments of colonocyte

rules = {'[lu]',  '[luLI]';
            '[luC]', '[luLI]';
           '[bpC]', '[e]'};

ColonModelMale_mergedCompsModel = mergeCompartmentsLIEC(ColonModelMale, rules);
ColonModelFemale_mergedCompsModel = mergeCompartmentsLIEC(ColonModelFemale, rules);

ColonModelMaleupd_mergedCompsModel = mergeCompartmentsLIEC(ColonModelMaleupd, rules);
ColonModelFemaleupd_mergedCompsModel = mergeCompartmentsLIEC(ColonModelFemaleupd, rules);

%Remaining compartments:
%    {'bc'  }
%    {'bp'  }
%    {'c'   }
%    {'d'   }
%    {'e'   }
%    {'g'   }
%    {'i'   }
%    {'l'   }
%    {'luLI'}
%    {'m'   }
%    {'n'   }
%    {'r'   }
%    {'x'   }

colonModelMale_withDietFe_mergedCompsModel = mergeCompartmentsLIEC(colonModelMale_withDietFe, rules);
colonModelFemale_withDietFe_mergedCompsModel = mergeCompartmentsLIEC(colonModelFemale_withDietFe, rules);

colonModelMale_withDietFeupd_mergedCompsModel = mergeCompartmentsLIEC(colonModelMale_withDietFeupd, rules);
colonModelFemale_withDietFeupd_mergedCompsModel = mergeCompartmentsLIEC(colonModelFemale_withDietFeupd, rules);

LIECmale = colonModelMale_withDietFeupd_mergedCompsModel;
LIECfemale = colonModelFemale_withDietFeupd_mergedCompsModel;

save(fullfile(modelsPath, 'LIECmale'),'LIECmale')
save(fullfile(modelsPath, 'LIECfemale'),'LIECfemale')

%Remaining compartments:
%    {'bc'  }
%    {'bp'  }
%    {'c'   }
%    {'d'   }
%    {'e'   }
%    {'g'   }
%    {'i'   }
%    {'l'   }
%    {'luLI'}
%    {'m'   }
%    {'n'   }
%    {'r'   }
%    {'x'   }

%% Test the models

objective   = 'biomass_maintenance';

% Define the name of the current model set (Change this for each run!)
setName = 'part04_mergeCompartmentsLIEC';
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmm')); 

% Resulting folder
resultsPath = fullfile(resultsPath, sprintf('%s_%s', setName, timestamp));

if ~exist(resultsPath, 'dir')
    mkdir(resultsPath);
end

% Pack Models into a Structure
% Load all your workspace models into a single structure.
% The field names used here will be saved in the master '00_INDEX' Excel sheet.
myModels = struct();

fprintf('Packing models into pipeline structure...\n');

% Base Models
myModels.ColonModelFemale = ColonModelFemale;
myModels.ColonModelMale = ColonModelMale;
myModels.ColonModelFemaleupd = ColonModelFemaleupd;
myModels.ColonModelMaleupd = ColonModelMaleupd;

% Merged Models
myModels.ColonModelFemale_Merged       = ColonModelFemale_mergedCompsModel;
myModels.colonModelMale_Merged         = ColonModelMale_mergedCompsModel;
myModels.colonModelFemaleupd_Merged    = ColonModelFemaleupd_mergedCompsModel;
myModels.colonModelMaleupd_Merged      = ColonModelMaleupd_mergedCompsModel;

% Diet Fe Models
myModels.colonModelFemale_DietFe       = colonModelFemale_withDietFe;
myModels.colonModelMale_DietFe         = colonModelMale_withDietFe;
myModels.colonModelFemale_DietFeMerged = colonModelFemale_withDietFe_mergedCompsModel;
myModels.colonModelMale_DietFeMerged   = colonModelMale_withDietFe_mergedCompsModel;

% Diet Fe Updated Models
myModels.colonFemale_DietFeUpd         = colonModelFemale_withDietFeupd;
myModels.colonMale_DietFeUpd           = colonModelMale_withDietFeupd;
myModels.colonFemale_DietFeUpdMerged   = colonModelFemale_withDietFeupd_mergedCompsModel;
myModels.colonMale_DietFeUpdMerged     = colonModelMale_withDietFeupd_mergedCompsModel;

%% 3. Run the Analysis Pipeline
% Call the function (defined below) to process everything and export to Excel
[PipelineResults, UpdatedModels] = analyzeMetabolicModels(myModels, objective, resultsPath, setName);

fprintf('\nPipeline completely finished! Check your results folder.\n');

%% 4. Run transport and exchange analysis

modelNames = fieldnames(myModels);

% Initialize the specific output structures you requested
transportSummary = struct();
transportMetsList = struct();
allRxns = struct();


for i = 1:length(modelNames)
    currentName = modelNames{i};
    currentModel = myModels.(currentName);
    
    fprintf('\nProcessing model %d of %d: %s\n', i, length(modelNames), currentName);
    
    % Run the function
    [transportSummary.(currentName), transportMetsList.(currentName), allRxns.(currentName)] = analyzeTransportsAndExchanges(currentModel);
    
end

fprintf('\nAll models processed! Check your workspace for the new variables.\n');


% Define the filenames for your three master Excel files
% (Make sure 'resultsPath' is defined in your workspace)
summaryFile = fullfile(resultsPath, 'Transport_Summaries.xlsx');
rxnsFile    = fullfile(resultsPath, 'All_Reactions.xlsx');
metsFile    = fullfile(resultsPath, 'Transport_Metabolites.xlsx');

% Get all the model names from your structure
modelNames = fieldnames(transportSummary);
numModels  = length(modelNames);

% Create a mapping array for the Master Index Sheet
sheetMapping = cell(numModels, 2);

fprintf('\nStarting Excel Export with Master Index...\n');

for i = 1:numModels
    currentModel = modelNames{i};
    
    % Create a perfectly safe, uniform sheet name (e.g., 'Model_001')
    safeSheetName = sprintf('Model_%03d', i);
    
    % Save to our mapping array
    sheetMapping{i, 1} = safeSheetName;
    sheetMapping{i, 2} = currentModel;
    
    fprintf('Exporting %s as %s...\n', currentModel, safeSheetName);
    
    % =========================================================
    % 1. File 1: Save the Transport Summary Table
    % =========================================================
    if ~isempty(transportSummary.(currentModel))
        writetable(transportSummary.(currentModel), summaryFile, 'Sheet', safeSheetName);
    end
    
    % =========================================================
    % 2. File 2: Save the All Reactions Table
    % =========================================================
    if ~isempty(allRxns.(currentModel))
        writetable(allRxns.(currentModel), rxnsFile, 'Sheet', safeSheetName);
    end
    
    % =========================================================
    % 3. File 3: Save the Metabolites List 
    % (Flattens the nested struct into a 2-column Excel table)
    % =========================================================
    metsStruct = transportMetsList.(currentModel);
    catNames = fieldnames(metsStruct);
    
    CategoryCol = {};
    MetaboliteCol = {};
    
    for k = 1:length(catNames)
        currentCat = catNames{k};
        mets = metsStruct.(currentCat);
        
        % Duplicate the category name to match the number of metabolites
        CategoryCol = [CategoryCol; repmat({currentCat}, length(mets), 1)];
        MetaboliteCol = [MetaboliteCol; mets];
    end
    
    % Convert to a table and write to Excel
    if ~isempty(CategoryCol)
        metsTable = table(CategoryCol, MetaboliteCol, 'VariableNames', {'Category', 'Metabolite'});
        writetable(metsTable, metsFile, 'Sheet', safeSheetName);
    end
end

% =========================================================
% 4. Write the Master Index Sheet to all files
% =========================================================
fprintf('\nWriting 00_INDEX master mapping sheets...\n');
T_Index = cell2table(sheetMapping, 'VariableNames', {'Excel_Sheet_Name', 'Full_MATLAB_Model_Name'});

% Suppress warnings about adding new sheets to Excel
warning('off', 'MATLAB:xlswrite:AddSheet');

writetable(T_Index, summaryFile, 'Sheet', '00_INDEX');
writetable(T_Index, rxnsFile, 'Sheet', '00_INDEX');
writetable(T_Index, metsFile, 'Sheet', '00_INDEX');

warning('on', 'MATLAB:xlswrite:AddSheet');

fprintf('✅ Excel export completed successfully!\n');