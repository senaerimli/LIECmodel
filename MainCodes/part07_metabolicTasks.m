%%metabolic tasks and core rxns
%by Sena Erimli

clear all
close all 
clc

initCobraToolbox(false)
addpath(genpath('.'));

% 1. DYNAMIC PATHS (Works on any computer)
scriptDir = fileparts(mfilename('fullpath')); 
basePath  = fullfile(scriptDir, 'LIEC'); 

global modelsPath resultsPath dataPath setName tasksPath
modelsPath  = fullfile(basePath, 'Models');
resultsPath = fullfile(basePath, 'Results');
dataPath    = fullfile(basePath, 'Data');
tasksPath   = fullfile(basePath, 'Tasks');
setName     = 'part07_metabolicTasks';
addpath(modelsPath, resultsPath, dataPath, tasksPath);

% Change solvers
LPSolver = 'gurobi';
QPSolver = 'gurobi';
changeCobraSolver(LPSolver, 'LP', 0, 1);
changeCobraSolver(QPSolver, 'QP', 0, 1);

% --- Models ---
fprintf('Loading draft reconstruction...\n');
%LIECfemale = readCbModel('LIECfemale.mat');
%LIECmale = readCbModel('LIECmale.mat');
mergedLIECfemale = readCbModel('mergedLIECfemale.mat');
mergedLIECmale = readCbModel('mergedLIECmale.mat');
%iColon = readCbModel('iColon.mat');

myModels = struct();

fprintf('Packing models into pipeline structure...\n');

% Base Models
%myModels.LIECfemale = LIECfemale;
%myModels.LIECmale = LIECmale;
myModels.mergedLIECfemale = mergedLIECfemale;
myModels.mergedLIECmale = mergedLIECmale;
%myModels.iColon = iColon;

% 1. Extract the names of all models stored in the structure
modelNames = fieldnames(myModels);

% Create an empty structure to store the results for each model
allTaskResults = struct();
sol = struct();


% 2. Loop through each model
for i = 1:length(modelNames)
    

    % Get the current model's name (e.g., 'mergedLIECfemale')
    modelID = modelNames{i};
    
    % Extract the actual COBRA model object using dynamic field referencing
    model = myModels.(modelID);
    
    fprintf('\n======================================================\n');
    fprintf('Evaluating Metabolic Tasks for: %s\n', modelID);
    fprintf('======================================================\n');
    
    metModel.(modelID) = runMetabolicTasks(model)

    %metModel.(modelID) = runMetabolicTasks(model, true)  %if core rxn
    %extention will be performed
end

for i = 1:length(modelNames)
    % Get the current model's name (e.g., 'mergedLIECfemale')
    modelID = modelNames{i};
    
    % Extract the actual COBRA model object using dynamic field referencing
    model = metModel.(modelID);

    objective   = 'biomass_maintenance';
    currentModel = changeObjective(model, objective);
    sol.(modelID) = optimizeCbModel(currentModel);
    disp(sol.(modelID))
end


LIECfemale = metModel.LIECfemale;
LIECmale = metModel.LIECmale;
mergedLIECfemale = metModel.mergedLIECfemale;
mergedLIECmale = metModel.mergedLIECmale;
iColon = metModel.iColon;

save(fullfile(modelsPath, 'LIECfemale'),'LIECfemale')
save(fullfile(modelsPath, 'LIECmale'),'LIECmale')
save(fullfile(modelsPath, 'mergedLIECfemale'),'mergedLIECfemale')
save(fullfile(modelsPath, 'mergedLIECmale'),'mergedLIECmale')
save(fullfile(modelsPath, 'iColon'),'iColon')

fprintf('\nPipeline complete.\n');