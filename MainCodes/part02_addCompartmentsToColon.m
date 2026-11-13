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

% Calling fileparts twice steps up from part02_addCompartmentsToColon.m -> MainCodes -> LIEC
basePath = fileparts(fileparts(scriptPath));

global modelsPath resultsPath 
modelsPath  = fullfile(basePath, 'Models');
resultsPath = fullfile(basePath, 'Results');

% Create the directories in LIEC if they don't exist
if ~exist(modelsPath, 'dir')
    mkdir(modelsPath);
end
if ~exist(resultsPath, 'dir')
    mkdir(resultsPath);
end

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

%%
%%DEFINE MODELS
male_unconstrained = loadPSCMfile('Harvey');
female_unconstrained = loadPSCMfile('Harvetta');

load OrganAtlas_Harvetta.mat;
load OrganAtlas_Harvey.mat;

ColonModelFemale = OrganCompendium_female.Colon.modelAllComp;
ColonModelMale = OrganCompendium_male.Colon.modelAllComp;

%load ColonModelFemale.mat;
%load ColonModelMale.mat;

%% ADD DIET COMPARTMENT
% Define your compartment tag (make sure it matches your model's format, e.g., '[c]' or '_c')
targetCompartment = '[d]';

% Find the reactions
% rxnList will contain the IDs, and rxnFormulaList will contain the equations
[d_rxnList_male] = findRxnFromCompartment(male_unconstrained, targetCompartment);
[d_rxnList_female] = findRxnFromCompartment(female_unconstrained, targetCompartment);

% Create a logical index that is 'true' wherever "Diet_" is found in the first column

hasGall_male = contains(d_rxnList_male(:, 2), 'Gall_');
hasGall_female = contains(d_rxnList_female(:, 2), 'Gall_');

% 2. Keep only the rows where hasGall is FALSE (using the ~ operator)
d_rxnList_male = d_rxnList_male(~hasGall_male, :);
d_rxnList_female = d_rxnList_female(~hasGall_female, :);

isDiet_male = contains(d_rxnList_male(:, 1), 'Diet_');
isDiet_female = contains(d_rxnList_female(:, 1), 'Diet_');

dietList_male = d_rxnList_male(isDiet_male, :);
dietList_female = d_rxnList_female(isDiet_female, :);

dietTrList_male = d_rxnList_male(~isDiet_male, :);
dietTrList_female = d_rxnList_female(~isDiet_female, :);


% Replace '[lu]' with '[luLI]' across the entire cell array instantly
updated_formulas_male = strrep(dietTrList_male(:, 2), '[lu]', '[luLI]');
updated_formulas_female = strrep(dietTrList_female(:, 2), '[lu]', '[luLI]');

% create submodel of Diet 
rxnsToRemove_male = setdiff(male_unconstrained.rxns, dietList_male);
dietSubModel_male = removeRxns(male_unconstrained, rxnsToRemove_male);

rxnsToRemove_female = setdiff(female_unconstrained.rxns, dietList_female);
dietSubModel_female = removeRxns(female_unconstrained, rxnsToRemove_female);

clearvars rxnsToRemove_male rxnsToRemove_female

% create submodel of diet transport reactions
rxnsToRemove_male = setdiff(male_unconstrained.rxns, dietTrList_male);
dietTrSubModel_male = removeRxns(male_unconstrained, rxnsToRemove_male);

rxnsToRemove_female = setdiff(female_unconstrained.rxns, dietTrList_female);
dietTrSubModel_female = removeRxns(female_unconstrained, rxnsToRemove_female);

[dietTrSubModel_male] = mergeCompartmentsLIEC(dietTrSubModel_male, {'[lu]',  '[luLI]'});
[dietTrSubModel_female] = mergeCompartmentsLIEC(dietTrSubModel_female, {'[lu]',  '[luLI]'});

% 1. Get the exact number of Fecal reactions you are adding
%numMets_male = length(dietSubModel_male.mets);
%numMets_female = length(dietSubModel_female.mets);

% 2. Create arrays of the exact same length for the bounds
%lb_male = -1000 * ones(numMets_male, 1);       % An array of all -1000s
%ub_male = zeros(numMets_male, 1);              % An array of all 0s
%lb_female = -1000 * ones(numMets_female, 1);       % An array of all -1000s
%ub_female = zeros(numMets_female, 1);              % An array of all 0s

% add Diet exchange reactions to the model -- make lb and ub == 0
%colonModelMale_withDiet = addExchangeRxn (ColonModelMale, dietSubModel_male.mets, lb_male, ub_male);
%colonModelFemale_withDiet = addExchangeRxn (ColonModelFemale, dietSubModel_female.mets, lb_female, ub_male);

colonModelMale_withDiet = addExchangeRxn (ColonModelMale, dietSubModel_male.mets, dietSubModel_male.lb, dietSubModel_male.ub);
colonModelFemale_withDiet = addExchangeRxn (ColonModelFemale, dietSubModel_female.mets, dietSubModel_female.lb, dietSubModel_female.ub);

% add diet transport rxns to the model
colonModelMale_withTrDiet = colonModelMale_withDiet;
colonModelFemale_withTrDiet = colonModelFemale_withDiet;

[compMets_male] = findMetFromCompartment(dietTrSubModel_male,'[luLI]');
[tf,dups] = checkIDsForTypeExist(colonModelMale_withDiet,compMets_male,'mets');   
clean_compMets_male = compMets_male(~ismember(compMets_male, dups));

clearvars tf dups

[compMets_female] = findMetFromCompartment(dietTrSubModel_female,'[luLI]');
[tf,dups] = checkIDsForTypeExist(colonModelFemale_withDiet,compMets_female,'mets');   
clean_compMets_female = compMets_female(~ismember(compMets_female, dups));

colonModelMale_withTrDiet = addMultipleMetabolites(colonModelMale_withTrDiet, clean_compMets_male);
colonModelFemale_withTrDiet = addMultipleMetabolites(colonModelFemale_withTrDiet, clean_compMets_female);

colonModelMale_withTrDietnew = colonModelMale_withTrDiet;
colonModelFemale_withTrDietnew = colonModelFemale_withTrDiet; 

for i = 1:length(dietTrSubModel_male.rxns)
    colonModelMale_withTrDietnew = addReaction(colonModelMale_withTrDietnew, ...
        dietTrSubModel_male.rxns{i}, ...
        'reactionFormula', updated_formulas_male{i}, ...
        'reactionName', dietTrSubModel_male.rxnNames{i}, ...
        'lowerBound', dietTrSubModel_male.lb(i), ...
        'upperBound', dietTrSubModel_male.ub(i), ...
        'subSystem', 'Transport Reaction'); 
end

for i = 1:length(dietTrSubModel_female.rxns)
    colonModelFemale_withTrDietnew = addReaction(colonModelFemale_withTrDietnew, ...
        dietTrSubModel_female.rxns{i}, ...
        'reactionFormula', updated_formulas_female{i}, ...
        'reactionName', dietTrSubModel_female.rxnNames{i}, ...
        'lowerBound', dietTrSubModel_female.lb(i), ...
        'upperBound', dietTrSubModel_female.ub(i), ...
        'subSystem', 'Transport Reaction'); 
end


%save(fullfile(modelsPath, 'male_ColonModelDiet_noConst.mat'), 'colonModelMale_withTrDietnew');
%save(fullfile(modelsPath, 'female_ColonModelDiet_noConst.mat'), 'colonModelFemale_withTrDietnew');
save(fullfile(modelsPath, 'male_ColonModelDiet.mat'), 'colonModelMale_withTrDietnew');
save(fullfile(modelsPath, 'female_ColonModelDiet.mat'), 'colonModelFemale_withTrDietnew');

%% NOW IT'S TIME FOR FECAL COMPARTMENT--- Not as new compartment 

excretionList_male = dietList_male;
excretionList_female = dietList_female;

excretionList_male(:, 1) = strrep(excretionList_male(:, 1), '[d]', '[luLI]');
excretionList_female(:, 1) = strrep(excretionList_female(:, 1), '[d]', '[luLI]');

excretionList_male(:, 2) = strrep(excretionList_male(:, 2), '[d]', '[luLI]');
excretionList_female(:, 2) = strrep(excretionList_female(:, 2), '[d]', '[luLI]');
excretionList_male(:, 2) = strrep(excretionList_male(:, 2), '->', '<=>');
excretionList_female(:, 2) = strrep(excretionList_female(:, 2), '->', '<=>');


clearvars compMets_female compMets_male clean_compMets_female clean_compMets_male tf dups

[compMets_male] = findMetFromCompartment(male_unconstrained,'[fe]');
compMets_male(:, 1) = strrep(compMets_male(:, 1), '[fe]', '[luLI]');

clearvars hasGall_female hasGall_male
hasGall_male = contains(compMets_male(:, 1), 'Gall_');

% 2. Keep only the rows where hasGall is FALSE (using the ~ operator)
compMets_male = compMets_male(~hasGall_male, :);

[tf,dups] = checkIDsForTypeExist(colonModelMale_withTrDietnew,compMets_male,'mets');   
clean_compMets_male = compMets_male(~ismember(compMets_male, dups));

clearvars tf dups

[compMets_female] = findMetFromCompartment(female_unconstrained,'[d]');
compMets_female(:, 1) = strrep(compMets_female(:, 1), '[d]', '[luLI]');

hasGall_female = contains(compMets_female(:, 1), 'Gall_');

% 2. Keep only the rows where hasGall is FALSE (using the ~ operator)
compMets_female = compMets_female(~hasGall_female, :);

[tf,dups] = checkIDsForTypeExist(colonModelMale_withTrDietnew,compMets_female,'mets');   
clean_compMets_female = compMets_female(~ismember(compMets_female, dups));

% 1. Get the exact number of Fecal reactions you are adding
numMets_male = size(compMets_male, 1);
numMets_female = size(compMets_female, 1);

% 2. Create arrays of the exact same length for the bounds
lb_male = zeros(numMets_male, 1);          % An array of all 0s
ub_male = 1000 * ones(numMets_male, 1);    % An array of all 1000s
lb_female = zeros(numMets_female, 1);          % An array of all 0s
ub_female = 1000 * ones(numMets_female, 1);    % An array of all 1000s

% add Fecal exchange reactions to the model -- make lb and ub == 0
colonModelMale_withDietFe = addExchangeRxn (colonModelMale_withTrDietnew, compMets_male, lb_male, ub_male);
colonModelFemale_withDietFe = addExchangeRxn (colonModelFemale_withTrDietnew, compMets_female, lb_female, ub_female);

%save(fullfile(modelsPath, 'male_ColonModelDietFe_noConst'),'colonModelMale_withDietFe')
%save(fullfile(modelsPath, 'female_ColonModelDietFe_noConst'),'colonModelFemale_withDietFe')
save(fullfile(modelsPath, 'male_ColonModelDietFe'),'colonModelMale_withDietFe')
save(fullfile(modelsPath, 'female_ColonModelDietFe'),'colonModelFemale_withDietFe')



%% ADD FECAL COMPARTMENT [fe] for every [d]

% 1. Replace '[d]' with '[fe]' at the exact end of all metabolite IDs
feTrSubModel_female = dietTrSubModel_female;
feTrSubModel_female.mets = regexprep(dietTrSubModel_female.mets, '\[d\]$', '[fe]');
feTrSubModel_female.rxns = regexprep(feTrSubModel_female.rxns, '\[d\]$', '[fe]');

% Update the model's compartment list if the 'comps' field exists
if isfield(feTrSubModel_female, 'comps')
    feTrSubModel_female.comps = regexprep(feTrSubModel_female.comps, '^d$', 'fe');
end

% Verify the changes
fprintf('Number of [fe] metabolites: %d\n', sum(endsWith(feTrSubModel_female.mets, '[fe]')));
fprintf('Number of [d] metabolites remaining: %d\n', sum(endsWith(feTrSubModel_female.mets, '[d]')));

% ==========================================
% 2. FIX: Properly chain the string replacements
% ==========================================
% Notice the fixed '[fe]' bracket and chaining the variable name
updated_formulas_female = strrep(dietTrList_female(:,2), '[d]', '[fe]');
updated_formulas_female = strrep(updated_formulas_female, '[lu]', '[luLI]');


% Initialize the new model
colonModelFemale_Fenew = colonModelFemale_withTrDietnew;

% ==========================================
% 3. FIX: Add Exchanges ONLY for the [fe] pool (With Array Bounds)
% ==========================================
% Find only the metabolites that end in [fe]
fe_mets_only = feTrSubModel_female.mets(endsWith(feTrSubModel_female.mets, '[fe]'));

% Generate boundary arrays of the exact same length (1354 elements)
num_fe_mets = length(fe_mets_only);
lb_array = zeros(num_fe_mets, 1);         % Array of 0s
ub_array = repmat(1000, num_fe_mets, 1);  % Array of 1000s


% Add the exchange reactions using the arrays
colonModelFemale_Fenew = addExchangeRxn(colonModelFemale_Fenew, fe_mets_only, lb_array, ub_array);

num_fe_mets = length(feTrSubModel_female.rxns);
lb_array = repmat(-1000, num_fe_mets, 1);  % Array of -1000s
ub_array = zeros(num_fe_mets, 1);         % Array of 0s

% ==========================================
% 4. Add the Transport Reactions
% ==========================================
for i = 1:length(feTrSubModel_female.rxns)
    colonModelFemale_Fenew = addReaction(colonModelFemale_Fenew, ...
        feTrSubModel_female.rxns{i}, ...
        'reactionFormula', updated_formulas_female{i}, ...
        'reactionName', feTrSubModel_female.rxnNames{i}, ...
        'lowerBound', lb_array(i), ...
        'upperBound', ub_array(i), ...
        'subSystem', 'Transport Reaction'); 
end

fprintf('Successfully integrated %d transport reactions and their fecal exchanges!\n', length(feTrSubModel_female.rxns));



%-----------
% Male model
%-----------

% 1. Replace '[d]' with '[fe]' at the exact end of all metabolite IDs
feTrSubModel_male = dietTrSubModel_male;
feTrSubModel_male.mets = regexprep(dietTrSubModel_male.mets, '\[d\]$', '[fe]');
feTrSubModel_male.rxns = regexprep(feTrSubModel_male.rxns, '\[d\]$', '[fe]');

% Update the model's compartment list if the 'comps' field exists
if isfield(feTrSubModel_male, 'comps')
    feTrSubModel_male.comps = regexprep(feTrSubModel_male.comps, '^d$', 'fe');
end

% Verify the changes
fprintf('Number of [fe] metabolites: %d\n', sum(endsWith(feTrSubModel_male.mets, '[fe]')));
fprintf('Number of [d] metabolites remaining: %d\n', sum(endsWith(feTrSubModel_male.mets, '[d]')));

% ==========================================
% 2. FIX: Properly chain the string replacements
% ==========================================
% Notice the fixed '[fe]' bracket and chaining the variable name
updated_formulas_male = strrep(dietTrList_male(:,2), '[d]', '[fe]');
updated_formulas_male = strrep(updated_formulas_male, '[lu]', '[luLI]');

% Initialize the new model
colonModelmale_Fenew = colonModelMale_withTrDietnew;

% ==========================================
% 3. FIX: Add Exchanges ONLY for the [fe] pool (With Array Bounds)
% ==========================================
% Find only the metabolites that end in [fe]
fe_mets_only = feTrSubModel_male.mets(endsWith(feTrSubModel_male.mets, '[fe]'));

% Generate boundary arrays of the exact same length (1354 elements)
num_fe_mets = length(fe_mets_only);
lb_array = zeros(num_fe_mets, 1);         % Array of 0s
ub_array = repmat(1000, num_fe_mets, 1);  % Array of 1000s

% Add the exchange reactions using the arrays
colonModelmale_Fenew = addExchangeRxn(colonModelmale_Fenew, fe_mets_only, lb_array, ub_array);


num_fe_mets = length(feTrSubModel_male.rxns);
lb_array = repmat(-1000, num_fe_mets, 1);  % Array of -1000s
ub_array = zeros(num_fe_mets, 1);         % Array of 0s

% ==========================================
% 4. Add the Transport Reactions
% ==========================================
for i = 1:length(feTrSubModel_male.rxns)
    colonModelmale_Fenew = addReaction(colonModelmale_Fenew, ...
        feTrSubModel_male.rxns{i}, ...
        'reactionFormula', updated_formulas_male{i}, ...
        'reactionName', feTrSubModel_male.rxnNames{i}, ...
        'lowerBound', lb_array(i), ...
        'upperBound', ub_array(i), ...
        'subSystem', 'Transport Reaction'); 
end

fprintf('Successfully integrated %d transport reactions and their fecal exchanges!\n', length(feTrSubModel_male.rxns));

save(fullfile(modelsPath, 'colonModelmale_Fenew'),'colonModelmale_Fenew')
save(fullfile(modelsPath, 'colonModelfemale_Fenew'),'colonModelFemale_Fenew')
