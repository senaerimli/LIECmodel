function [metModel, resultVac, missedVac, coreRxnsVac] = runMetabolicTasks(model, doCoreExtraction)
%% runMetabolicTasks.m
% Pipeline runner for strict vacuum metabolic task testing and core reaction extraction.

if nargin < 2
    doCoreExtraction = false;
end

%% 1. Load reference model and define solver options
fprintf('Loading reference model (Recon3D)...\n');
generalModel = readCbModel('Recon3DModel_301.mat');

TASKS_FILE = 'ColonMetabolicTasks_final.xlsx'; 

opts.solver          = 'gurobi';    
opts.tol             = 1e-6;        
opts.verbose         = true;        
opts.extractCoreRxns = doCoreExtraction; 

%% 2. Output directory 
global setName resultsPath
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmm')); 
OUTPUT_DIR = fullfile(resultsPath, sprintf('%s_%s_%s', setName, model.modelID, timestamp));

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% 3. Run Metabolic Tasks in Strict Vacuum (Topology & Core Diagnostics)
fprintf('\n======================================================\n');
fprintf(' METABOLIC TASK TESTING (Strict Vacuum - %s)\n', model.modelID);
fprintf('======================================================\n');

% Step 1: Enforce internal colonocyte tissue-specific restrictions FIRST
model_Colon = updateColonSpecBoundaries(model);

% Step 2: Seal external boundaries to a strict vacuum SECOND
model_Vac = vacuumTheModel(model_Colon);

% Step 3: Evaluate all metabolic tasks
opts.outputFile = fullfile(OUTPUT_DIR, 'Results_Vacuum_Tasks.xlsx');
[resultVac, missedVac, coreRxnsVac] = testColonMetabolicTasks(model_Vac, generalModel, TASKS_FILE, opts);

%% 4. Diagnostic Summary (Positive Tasks vs. [SHOULD FAIL] Negative Controls)
isShouldFail = contains(upper(resultVac.TaskName), 'FAIL');

nPosTotal = sum(~isShouldFail);
nPosPass  = sum(resultVac.Passed(~isShouldFail));

nNegTotal = sum(isShouldFail);
nNegFail  = sum(~resultVac.Passed(isShouldFail));

fprintf('\n======================================================\n');
fprintf(' VACUUM TASK SUMMARY (%s)\n', model.modelID);
fprintf('======================================================\n');
fprintf('  Physiological Tasks Passed     : %d / %d (%.1f%%)\n', nPosPass, nPosTotal, 100 * nPosPass / max(nPosTotal, 1));
fprintf('  [SHOULD FAIL] Controls Failed  : %d / %d (%.1f%%)\n', nNegFail, nNegTotal, 100 * nNegFail / max(nNegTotal, 1));
fprintf('  Results saved to               : %s\n', opts.outputFile);

%% 5. Save Essential Core Reactions (if requested)
if doCoreExtraction
    coreFile = fullfile(OUTPUT_DIR, sprintf('CoreReactions_Vacuum_%s.mat', model.modelID));
    save(coreFile, 'coreRxnsVac');
    fprintf('  Essential Core Reactions saved : %s\n', coreFile);
end

% Return the colon-constrained model (with tissue boundaries applied, prior to vacuuming)
metModel = model_Colon;
end