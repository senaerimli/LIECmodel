%%WBM extraction
%by Sena Erimli
% this code is actually a snippet from Harvey/Harvetta by Thiele group. I
% have used their published code to identify boundary differences between
% the models with different parameters

%[OrganCompendium_male,TableCSources] = getOrgansFromHarvey(male);
%    OrganCompendium_male.sex = 'male';

 %[OrganCompendium_female,TableCSources] = getOrgansFromHarvey(female);
 %   OrganCompendium_female.sex = 'female';

%Female_b4_32=OrganCompendium_female.Colon.modelAllComp;
%Male_b4_32=OrganCompendium_male.Colon.modelAllComp;

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
%You can skip here if you have saved files
male = loadPSCMfile('Harvey');
female = loadPSCMfile('Harvetta');

male_unconstrained = male;
female_unconstrained = female;

%save(fullfile(modelsPath, 'male_unconstrained'),'male_unconstrained')
%save(fullfile(modelsPath, 'female_unconstrained'),'female_unconstrained')

% For the male:

sex = 'male';
standardPhysiolDefaultParameters;
IndividualParameters_male = IndividualParameters

% and for the female:

sex = 'female';
standardPhysiolDefaultParameters;
IndividualParameters_female = IndividualParameters

% apply HMDB metabolomic data based on the female individual parameters
female = physiologicalConstraintsHMDBbased(female,IndividualParameters_female);
male = physiologicalConstraintsHMDBbased(male,IndividualParameters_male);

male_physConstsHMDBbased = male;
female_physConstsHMDBbased = female;

%save(fullfile(modelsPath, 'male_physConstsHMDBbased'),'male_physConstsHMDBbased')
%save(fullfile(modelsPath, 'female_physConstsHMDBbased'),'female_physConstsHMDBbased')


%female
% set some more constraints
female = setSimulationConstraints(female);

female.lb(strmatch('BBB_KYNATE[CSF]upt',female.rxns)) = -1000000;
female.lb(strmatch('BBB_LKYNR[CSF]upt',female.rxns)) = -1000000;
female.lb(strmatch('BBB_TRP_L[CSF]upt',female.rxns)) = -1000000;
female.ub(strmatch('Brain_EX_glc_D(',female.rxns)) = -100; % currently -400 rendering many of the models to be infeasible in germfree state

minInf = -1000000;
maxInf = 1000000;
minConstraints = length(intersect(find(female.lb>minInf),find(female.lb)));
maxConstraints =length(intersect(find(female.ub<maxInf),find(female.ub)));
PercentageConstraintRxns_female = (minConstraints + maxConstraints)*100/length(female.ub)
% male
% set some more constraints
male = setSimulationConstraints(male);
male.lb(strmatch('BBB_KYNATE[CSF]upt',male.rxns)) = -1000000; %constrained uptake
male.lb(strmatch('BBB_LKYNR[CSF]upt',male.rxns)) = -1000000; %constrained uptake
male.lb(strmatch('BBB_TRP_L[CSF]upt',male.rxns)) = -1000000; %constrained uptake
male.ub(strmatch('Brain_EX_glc_D(',male.rxns)) = -100; % currently -400 rendering many of the models to be infeasible in germfree state

% count how many reactions have a non-infinity bounds
minInf = -1000000;
maxInf = 1000000;
minConstraints = length(intersect(find(male.lb>minInf),find(male.lb)));
maxConstraints =length(intersect(find(male.ub<maxInf),find(male.ub)));
PercentageConstraintRxns_male = (minConstraints + maxConstraints)*100/length(male.ub)

male_ConstsWoDiet = male;
female_ConstsWoDiet = female;

%save(fullfile(modelsPath, 'male_ConstsWoDiet'),'male_ConstsWoDiet')
%save(fullfile(modelsPath, 'female_ConstsWoDiet'),'female_ConstsWoDiet')

%add diet constraints
EUAverageDietNew

female = setDietConstraints(female, Diet);
male = setDietConstraints(male, Diet);

male_Consts = male; %Harvey_1_03c
female_Consts = female; %Harvetta_1_03c

%save(fullfile(modelsPath, 'male_Consts'),'male_Consts')
%save(fullfile(modelsPath, 'female_Consts'),'female_Consts')

%% --- Calculate Constrained Reactions for All Snapshots (Original Logic) ---

% Group the models into arrays for easy looping
models_male = {male_unconstrained, male_physConstsHMDBbased, male_ConstsWoDiet, male_Consts};
models_female = {female_unconstrained, female_physConstsHMDBbased, female_ConstsWoDiet, female_Consts};
model_names = {'unconstrained', 'physConstsHMDBbased', 'ConstsWoDiet', 'Consts'};

minInf = -1000000;
maxInf = 1000000;

fprintf('\n--- Values for Table 2 ---\n');
for i = 1:4
    % Male calculation
    m = models_male{i};
    minConstraints_m = length(intersect(find(m.lb>minInf),find(m.lb)));
    maxConstraints_m = length(intersect(find(m.ub<maxInf),find(m.ub)));
    pct_m = (minConstraints_m + maxConstraints_m)*100/length(m.ub);
    
    % Female calculation
    f = models_female{i};
    minConstraints_f = length(intersect(find(f.lb>minInf),find(f.lb)));
    maxConstraints_f = length(intersect(find(f.ub<maxInf),find(f.ub)));
    pct_f = (minConstraints_f + maxConstraints_f)*100/length(f.ub);
    
    % Print formatted output
    fprintf('\\texttt{*_%s} & [%.2f\\%% / %.2f\\%%] \\\\\n', model_names{i}, pct_m, pct_f);
end
fprintf('--------------------------\n');