%%Boundary comparison
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

%% Colon specific model and ub&lb comparison with other models

load male_unconstrained.mat;
load female_unconstrained.mat;

load OrganAtlas_Harvetta.mat;
load OrganAtlas_Harvey.mat;

load male_Consts.mat;
load female_Consts.mat;

ColonModelFemale=OrganCompendium_female.Colon.modelAllComp;
ColonModelMale=OrganCompendium_male.Colon.modelAllComp;

%save(fullfile(modelsPath, 'ColonModelMale'),'ColonModelMale')
%save(fullfile(modelsPath, 'ColonModelFemale'),'ColonModelFemale')

%compile part_02_addCompartmentsToColon.m
load male_ColonModelDietFe_noConst   %colonModelMale_withDietFe
load female_ColonModelDietFe_noConst %colonModelFemale_withDietFe

%comparing any of the model is enough
%first extract the related rxns. since rxn names are different, comparison
%will be performed via S matrix

compartments = {'[luLI]', '[fe]', '[bpC]', '[d]', '[luC]'};
prefixes = {'Colon_'};

subModel_male_unconstrained = extractCompartmentsFromModel(male_unconstrained, compartments, prefixes, true);
subModel_female_unconstrained = extractCompartmentsFromModel(female_unconstrained, compartments, prefixes, true);

%% SKIP THIS PART
% if you have obtained diet constraints' table

% 1. Find reactions ONLY in the Male model
rxnsOnlyMale = setdiff(ColonModelMale.rxns, ColonModelFemale.rxns);

% 2. Find reactions ONLY in the Female model
rxnsOnlyFemale = setdiff(ColonModelFemale.rxns, ColonModelMale.rxns);

% 3. Find reactions they both share (Intersection)
sharedRxns = intersect(ColonModelMale.rxns, ColonModelFemale.rxns);

% 1. Find the exact indices of the unique reactions in the Male model
[~, idxM] = ismember(rxnsOnlyMale, ColonModelMale.rxns);
[~, idxF] = ismember(rxnsOnlyFemale, ColonModelFemale.rxns);

% 2. Extract the required fields using those indices
rxns_M = ColonModelMale.rxns(idxM);
rxns_F = ColonModelFemale.rxns(idxF);
reactions_M = ColonModelMale.reactions(idxM);
reactions_F = ColonModelFemale.reactions(idxF);
subSystems_M = ColonModelMale.subSystems(idxM);
subSystems_F = ColonModelFemale.subSystems(idxF);

T_UniqueMale = table(rxns_M, reactions_M, subSystems_M, ...
    'VariableNames', {'rxns', 'reactions', 'subSystems'});

T_UniqueFemale = table(rxns_F, reactions_F, subSystems_F, ...
    'VariableNames', {'rxns', 'reactions', 'subSystems'});

%exportPath = fullfile(resultsPath, 'Part03_Isolated_Reactions_Details.xlsx');
%writetable(T_UniqueMale, exportPath, 'Sheet', 'Male_Unique');
%writetable(T_UniqueFemale, exportPath, 'Sheet', 'Female_Unique');


load male_physConstsHMDBbased.mat;
load female_physConstsHMDBbased.mat;

load male_ConstsWoDiet.mat;
load female_ConstsWoDiet.mat;

load colonModelmale_Fenew.mat;
load colonModelfemale_Fenew.mat;


subModel_male_Consts = extractCompartmentsFromModel(male_Consts, compartments, prefixes, true);
subModel_female_Consts = extractCompartmentsFromModel(female_Consts, compartments, prefixes, true);

%subModel_male_ConstsWoDiet = extractCompartmentsFromModel(male_ConstsWoDiet, compartments, prefixes, true);
%subModel_female_ConstsWoDiet = extractCompartmentsFromModel(female_ConstsWoDiet, compartments, prefixes, true);

%subModel_male_physConstsHMDBbased = extractCompartmentsFromModel(male_physConstsHMDBbased, compartments, prefixes, true);
%subModel_female_physConstsHMDBbased = extractCompartmentsFromModel(female_physConstsHMDBbased, compartments, prefixes, true);
%check the ub&lb of WBMs/ all identical
%Diff_unconstrained_physConstsHMDBbased_male= findBoundDiff(subModel_male_physConstsHMDBbased, subModel_male_unconstrained);
%Diff_unconstrained_physConstsHMDBbased_female= findBoundDiff(subModel_female_physConstsHMDBbased, subModel_female_unconstrained);
%all identical bounds
%Diff_unconstrained_ConstsWoDiet_male= findBoundDiff(subModel_male_ConstsWoDiet, subModel_male_unconstrained);
%Diff_unconstrained_ConstsWoDiet_female= findBoundDiff(subModel_female_ConstsWoDiet, subModel_female_unconstrained);
%all identical bounds

Diff_unconstrained_Consts_male= findBoundDiff(subModel_male_Consts, subModel_male_unconstrained);
Diff_unconstrained_Consts_female= findBoundDiff(subModel_female_Consts, subModel_female_unconstrained);
%36 rxns-Male
%36 rxns-Female

%36 diet exhange-Male
%36 diet exhange-Female
%T_diff_unConsts_male = compareExchanges(subModel_male_unconstrained, subModel_male_Consts, 'Unconst', 'Consts');
%T_diff_unConsts_female = compareExchanges(subModel_female_unconstrained, subModel_female_Consts, 'Unconst', 'Consts');
%36 diet exhange-Male
%36 diet exhange-Female
%same as Diff_unconstrained_Consts_male & Diff_unconstrained_Consts_female


% Find the common IDs
% 'commonIDs' will contain the list of identical reactions
%commonIDs_unConsts_male = intersect(Diff_unconstrained_Consts_male{:, 1}, T_diff_unConsts_male{:, 1});
%commonIDs_unConsts_female = intersect(Diff_unconstrained_Consts_female{:, 1}, T_diff_unConsts_female{:, 1});
% ALL IDENTICAL

diet_unconstrained_male.rxns = Diff_unconstrained_Consts_male{:, 2};
diet_Consts_male.rxns = Diff_unconstrained_Consts_male{:, 1};
diet_unconstrained_male.lb = Diff_unconstrained_Consts_male{:, 5};
diet_Consts_male.lb = Diff_unconstrained_Consts_male{:, 4};
diet_unconstrained_male.ub = Diff_unconstrained_Consts_male{:, 7};
diet_Consts_male.ub = Diff_unconstrained_Consts_male{:, 6};

diet_unconstrained_female.rxns = Diff_unconstrained_Consts_female{:, 2};
diet_Consts_female.rxns = Diff_unconstrained_Consts_female{:, 1};
diet_unconstrained_female.lb = Diff_unconstrained_Consts_female{:, 5};
diet_Consts_female.lb = Diff_unconstrained_Consts_female{:, 4};
diet_unconstrained_female.ub = Diff_unconstrained_Consts_female{:, 7};
diet_Consts_female.ub = Diff_unconstrained_Consts_female{:, 6};

% --- Create Male Diet Comparison Table ---
T_Diet = table(diet_Consts_male.rxns, diet_unconstrained_male.rxns, ...
    diet_Consts_male.lb, diet_unconstrained_male.lb, ...
    diet_Consts_male.ub, diet_unconstrained_male.ub, ...
    diet_Consts_female.rxns, diet_unconstrained_female.rxns, ...
    diet_Consts_female.lb, diet_unconstrained_female.lb, ...
    diet_Consts_female.ub, diet_unconstrained_female.ub, ...
    'VariableNames', {'Male_Const_Rxn', 'Male_Unconst_Rxn', 'Male_Const_LB',...
    'Male_Unconst_LB', 'Male_Const_UB', 'Male_Unconst_UB', 'Female_Const_Rxn',...
    'Female_Unconst_Rxn', 'Female_Const_LB', 'Female_Unconst_LB', ...
    'Female_Const_UB', 'Female_Unconst_UB'});


%writetable(T_Diet,  fullfile(resultsPath, 'Part03_Diet_Constraints.xlsx'));

%T_diff_unHMDB = compareExchanges(subModel_male_unconstrained, subModel_male_physConstsHMDBbased, 'Unconst', 'HMDB');
%nodiff
%T_diff_unConstsWoDiet = compareExchanges(subModel_male_unconstrained, subModel_male_ConstsWoDiet, 'Unconst', 'ConstWoDiet');
%nodiff



%%
% ----- CONTINUE FROM HERE -----

%Diff_ColonMale_unconstrained = findBoundDiff(ColonModelMale, subModel_male_unconstrained);
%Diff_ColonFemale_unconstrained = findBoundDiff(ColonModelFemale, subModel_female_unconstrained);
%704 rxns-Male
%700 rxns-Female

Diff_ColonMaleDietFe_unconstrained = findBoundDiff(colonModelMale_withDietFe, subModel_male_unconstrained);
Diff_ColonFemaleDietFe_unconstrained = findBoundDiff(colonModelFemale_withDietFe, subModel_female_unconstrained);
%2057 rxns-Male
%2054 rxns-Female

%Diff_ColonMale_DietFe = findBoundDiff(ColonModelMale, colonModelMale_withDietFe);
%Diff_ColonFemale_DietFe = findBoundDiff(ColonModelFemale, colonModelFemale_withDietFe);
% all identical

% Find the common IDs
% 'commonIDs' will contain the list of identical reactions
%commonIDs_ColonMale_DietFe = intersect(Diff_ColonMale_unconstrained{:, 1}, Diff_ColonMaleDietFe_unconstrained{:, 1});
%commonIDs_ColonFemale_DietFe = intersect(Diff_ColonFemale_unconstrained{:, 1}, Diff_ColonFemaleDietFe_unconstrained{:, 1});

% all rxns in Diff_ColonMale_unconstrained is in Diff_ColonMaleDietFe_unconstrained

%% --- EXTRACT [luC], [bc] and [luLI] BOUNDARIES ---
idsM_luCbC = Diff_ColonMaleDietFe_unconstrained{:, 1};
idsF_luCbC = Diff_ColonFemaleDietFe_unconstrained{:, 1};

maskM_luCbC = (contains(idsM_luCbC, '[luC]') | (contains(idsM_luCbC, '[bc]') & ~contains(idsM_luCbC, '_[bc]')));
maskF_luCbC = (contains(idsF_luCbC, '[luC]') | (contains(idsF_luCbC, '[bc]') & ~contains(idsF_luCbC, '_[bc]')));

T_filtered_male_luCbC = Diff_ColonMaleDietFe_unconstrained(maskM_luCbC, [1, 5, 7]);
T_filtered_male_luCbC.Properties.VariableNames = {'Reaction', 'maleM2_lb', 'maleM2_ub'};

T_filtered_female_luCbC = Diff_ColonFemaleDietFe_unconstrained(maskF_luCbC, [1, 5, 7]);
T_filtered_female_luCbC.Properties.VariableNames = {'Reaction', 'femaleM2_lb', 'femaleM2_ub'};

% --- EXTRACT [d] BOUNDARIES ---
rowMask_d_male = contains(Diff_ColonMaleDietFe_unconstrained{:, 1}, '[d]');
rowMask_d_female = contains(Diff_ColonFemaleDietFe_unconstrained{:, 1}, '[d]');

T_d_boundaries_male = Diff_ColonMaleDietFe_unconstrained(rowMask_d_male, [1, 5, 7]);
T_d_boundaries_female = Diff_ColonFemaleDietFe_unconstrained(rowMask_d_female, [1, 5, 7]);

T_d_boundaries_male.Properties.VariableNames = {'Reaction', 'maleM2_lb', 'maleM2_ub'};
T_d_boundaries_female.Properties.VariableNames = {'Reaction', 'femaleM2_lb', 'femaleM2_ub'};

%% --- Update Upper Bounds for ColonModel Model ---

% update [luC] and [bc]  boundaries
ColonModelMaleupd = updateModelBounds(ColonModelMale, T_filtered_male_luCbC, 'maleM2_ub');
ColonModelMaleupd = updateModelBounds(ColonModelMaleupd, T_filtered_male_luCbC, 'maleM2_lb');

ColonModelFemaleupd = updateModelBounds(ColonModelFemale, T_filtered_female_luCbC, 'femaleM2_ub');
ColonModelFemaleupd = updateModelBounds(ColonModelFemaleupd, T_filtered_female_luCbC, 'femaleM2_lb');

colonModelMale_Feupd = updateModelBounds(colonModelmale_Fenew, T_filtered_male_luCbC, 'maleM2_ub');
colonModelMale_Feupd = updateModelBounds(colonModelMale_Feupd, T_filtered_male_luCbC, 'maleM2_lb');

colonModelFemale_Feupd = updateModelBounds(colonModelFemale_Fenew, T_filtered_female_luCbC, 'femaleM2_ub');
colonModelFemale_Feupd = updateModelBounds(colonModelFemale_Feupd, T_filtered_female_luCbC, 'femaleM2_lb');


%save(fullfile(modelsPath, 'ColonModelMaleupd'),'ColonModelMaleupd')
%save(fullfile(modelsPath, 'ColonModelFemaleupd'),'ColonModelFemaleupd')

% update [d] boundary (ub 0, lb -1000000)
d_boundaries_male = T_d_boundaries_male;
d_boundaries_female = T_d_boundaries_female;

d_boundaries_male{:, 'maleM2_ub'} = 0;
d_boundaries_male{:, 'maleM2_lb'} = -1000000; 
d_boundaries_female{:, 'femaleM2_ub'} = 0;
d_boundaries_female{:, 'femaleM2_lb'} = -1000000;

colonModelMale_withDietFeupd = updateModelBounds(colonModelMale_withDietFe, d_boundaries_male, 'maleM2_ub');
colonModelMale_withDietFeupd = updateModelBounds(colonModelMale_withDietFeupd, d_boundaries_male, 'maleM2_lb');

colonModelFemale_withDietFeupd = updateModelBounds(colonModelFemale_withDietFe, d_boundaries_female, 'femaleM2_ub');
colonModelFemale_withDietFeupd = updateModelBounds(colonModelFemale_withDietFeupd, d_boundaries_female, 'femaleM2_lb');

colonModelMale_Feupd = updateModelBounds(colonModelmale_Fenew, d_boundaries_male, 'maleM2_ub');
colonModelMale_Feupd = updateModelBounds(colonModelMale_Feupd, d_boundaries_male, 'maleM2_lb');

colonModelFemale_Feupd = updateModelBounds(colonModelFemale_Fenew, d_boundaries_female, 'femaleM2_ub');
colonModelFemale_Feupd = updateModelBounds(colonModelFemale_Feupd, d_boundaries_female, 'femaleM2_lb');

colonModelMale_withDietFeupd = colonModelMale_Feupd;
colonModelFemale_withDietFeupd = colonModelFemale_Feupd;
%save(fullfile(modelsPath, 'colonModelMale_withDietFeupd'),'colonModelMale_withDietFeupd')
%save(fullfile(modelsPath, 'colonModelFemale_withDietFeupd'),'colonModelFemale_withDietFeupd')

%% Test the models

fprintf('\n--- Dead End Analysis ---\n');
[DeadIdx_male] = detectDeadEnds(ColonModelMale);
[newDeadIdx_male] = detectDeadEnds(ColonModelMaleupd);
fprintf('Total Dead Ends before constraints: %d\n', length(DeadIdx_male));
fprintf('Total Dead Ends after constraints: %d\n', length(newDeadIdx_male));

[DeadIdx_female] = detectDeadEnds(ColonModelFemale);
[newDeadIdx_female] = detectDeadEnds(ColonModelFemaleupd);
fprintf('Total Dead Ends before constraints: %d\n', length(DeadIdx_female));
fprintf('Total Dead Ends after constraints: %d\n', length(newDeadIdx_female));

[DeadIdx_df_male] = detectDeadEnds(colonModelMale_withDietFe);
[newDeadIdx_df_male] = detectDeadEnds(colonModelMale_withDietFeupd);
fprintf('Total Dead Ends before constraints: %d\n', length(DeadIdx_df_male));
fprintf('Total Dead Ends after constraints: %d\n', length(newDeadIdx_df_male));

[DeadIdx_df_female] = detectDeadEnds(colonModelFemale_withDietFe);
[newDeadIdx_df_female] = detectDeadEnds(colonModelFemale_withDietFeupd);

fprintf('Total Dead Ends before constraints: %d\n', length(DeadIdx_df_female));
fprintf('Total Dead Ends after constraints: %d\n', length(newDeadIdx_df_female));

fprintf('\n--- FBA Optimization ---\n');
sol = optimizeCbModel(ColonModelFemale)
sol = optimizeCbModel(ColonModelFemaleupd)

sol = optimizeCbModel(ColonModelMale)
sol = optimizeCbModel(ColonModelMaleupd)

sol = optimizeCbModel(colonModelFemale_withDietFe)
sol = optimizeCbModel(colonModelFemale_withDietFeupd)

sol = optimizeCbModel(colonModelMale_withDietFe)
sol = optimizeCbModel(colonModelMale_withDietFeupd)

%% DIRECT MODEL BOUNDARY EXTRACTION (DIET & EXCHANGE)
fprintf('\nExtracting and merging boundaries directly from models...\n');

% The Strict Universal Cleaner
cleanUniv = @(x) regexprep(x, '^(Diet_|Gall_|Excretion_|GI_|LI_|SI_|Tr_)', '');

% =========================================================================
% PART A: DIET BOUNDARIES TABLE
% =========================================================================
% Find Diet reactions (WBMs use 'Diet_', Colon models use '[d]')
idx_diet_unconst_m = contains(male_unconstrained.rxns, 'Diet_');
idx_diet_const_m   = contains(male_Consts.rxns, 'Diet_');
idx_diet_colon_m   = contains(colonModelMale_withDietFe.rxns, '[d]');

idx_diet_unconst_f = contains(female_unconstrained.rxns, 'Diet_');
idx_diet_const_f   = contains(female_Consts.rxns, 'Diet_');
idx_diet_colon_f   = contains(colonModelFemale_withDietFe.rxns, '[d]');

% Extract Tables using the helper function
T_D_Un_M  = extractAndClean(male_unconstrained, idx_diet_unconst_m, 'Unconst_Male', cleanUniv);
T_D_Con_M = extractAndClean(male_Consts, idx_diet_const_m, 'Const_Male', cleanUniv);
T_D_Col_M = extractAndClean(colonModelMale_withDietFe, idx_diet_colon_m, 'Colon_Male', cleanUniv);

T_D_Un_F  = extractAndClean(female_unconstrained, idx_diet_unconst_f, 'Unconst_Female', cleanUniv);
T_D_Con_F = extractAndClean(female_Consts, idx_diet_const_f, 'Const_Female', cleanUniv);
T_D_Col_F = extractAndClean(colonModelFemale_withDietFe, idx_diet_colon_f, 'Colon_Female', cleanUniv);

% Iteratively merge into a single Diet Master Table
T_Diet_Master = outerjoin(T_D_Un_M, T_D_Con_M, 'Keys', 'BaseMetabolite', 'MergeKeys', true);
T_Diet_Master = outerjoin(T_Diet_Master, T_D_Col_M, 'Keys', 'BaseMetabolite', 'MergeKeys', true);
T_Diet_Master = outerjoin(T_Diet_Master, T_D_Un_F, 'Keys', 'BaseMetabolite', 'MergeKeys', true);
T_Diet_Master = outerjoin(T_Diet_Master, T_D_Con_F, 'Keys', 'BaseMetabolite', 'MergeKeys', true);
T_Diet_Master = outerjoin(T_Diet_Master, T_D_Col_F, 'Keys', 'BaseMetabolite', 'MergeKeys', true);

%writetable(T_Diet_Master, fullfile(resultsPath, 'Part03_Direct_Diet_Comparison.xlsx'));
%fprintf('Saved Direct_Diet_Comparison.xlsx\n');

% =========================================================================
% PART B: EXCHANGE BOUNDARIES TABLE
% =========================================================================
% Use COBRA's native findExcRxns to isolate true exchange boundaries
idx_ex_unconst_m = findExcRxns(male_unconstrained);
idx_ex_const_m   = findExcRxns(male_Consts);
idx_ex_colon_m   = findExcRxns(colonModelMale_withDietFe);

idx_ex_unconst_f = findExcRxns(female_unconstrained);
idx_ex_const_f   = findExcRxns(female_Consts);
idx_ex_colon_f   = findExcRxns(colonModelFemale_withDietFe);

% Extract Tables
T_E_Un_M  = extractAndClean(male_unconstrained, idx_ex_unconst_m, 'Unconst_Male', cleanUniv);
T_E_Con_M = extractAndClean(male_Consts, idx_ex_const_m, 'Const_Male', cleanUniv);
T_E_Col_M = extractAndClean(colonModelMale_withDietFe, idx_ex_colon_m, 'Colon_Male', cleanUniv);

T_E_Un_F  = extractAndClean(female_unconstrained, idx_ex_unconst_f, 'Unconst_Female', cleanUniv);
T_E_Con_F = extractAndClean(female_Consts, idx_ex_const_f, 'Const_Female', cleanUniv);
T_E_Col_F = extractAndClean(colonModelFemale_withDietFe, idx_ex_colon_f, 'Colon_Female', cleanUniv);

% Iteratively merge into a single Exchange Master Table
T_Ex_Master = outerjoin(T_E_Un_M, T_E_Con_M, 'Keys', 'BaseMetabolite', 'MergeKeys', true);
T_Ex_Master = outerjoin(T_Ex_Master, T_E_Col_M, 'Keys', 'BaseMetabolite', 'MergeKeys', true);
T_Ex_Master = outerjoin(T_Ex_Master, T_E_Un_F, 'Keys', 'BaseMetabolite', 'MergeKeys', true);
T_Ex_Master = outerjoin(T_Ex_Master, T_E_Con_F, 'Keys', 'BaseMetabolite', 'MergeKeys', true);
T_Ex_Master = outerjoin(T_Ex_Master, T_E_Col_F, 'Keys', 'BaseMetabolite', 'MergeKeys', true);

%writetable(T_Ex_Master, fullfile(resultsPath, 'Part03_Direct_Exchange_Comparison.xlsx'));
%fprintf('Saved Direct_Exchange_Comparison.xlsx\n');

%% DEFINE TRANSPORT REACTIONS FOR CLARITY

%find transport rxns with upper and lower bounds different than 1000, -1000
%and 0
bc_transRxnsT_male_unconstrained = findTransportRxns(subModel_male_Consts, '[e]', '[bc]', false);
bc_transRxnsT_female_unconstrained = findTransportRxns(subModel_female_Consts, '[e]', '[bc]', false);
%321 rxns female (+3 rxns with 0 to 1000 bound)
%326 rxns male (+3 rxns with 0 to 1000 bound)

bp_transRxnsT_male_unconstrained = findTransportRxns(male_unconstrained, '[bpC]', '[bp]', false);
bp_transRxnsT_female_unconstrained = findTransportRxns(female_unconstrained, '[bpC]', '[bp]', false);
%412 rxns have 0 to 1000 bound

d_transRxnsT_male_unconstrained = findTransportRxns(male_unconstrained, '[d]', '[luLI]', false);
d_transRxnsT_female_unconstrained = findTransportRxns(female_unconstrained, '[d]', '[luLI]', false);
%1391 (female) and 1392 (male) rxns have 0 to 1000 bound

lu_transRxnsT_male_unconstrained = findTransportRxns(male_unconstrained, '[lu]', '[luSI]', false);
lu_transRxnsT_female_unconstrained = findTransportRxns(female_unconstrained, '[lu]', '[luSI]', false);
%1393 (female) and 1395 (male) rxns have 0 to 1000 bound

lu_LI_transRxnsT_male_unconstrained = findTransportRxns(male_unconstrained, '[lu]', '[luLI]', false);
lu_LI_transRxnsT_female_unconstrained = findTransportRxns(female_unconstrained, '[lu]', '[luLI]', false);
%no matching rxns

lu_C_transRxnsT_male_unconstrained = findTransportRxns(male_unconstrained, '[lu]', '[luC]', false);
lu_C_transRxnsT_female_unconstrained = findTransportRxns(female_unconstrained, '[lu]', '[luC]', false);
%no matching rxns

luSI_transRxnsT_male_unconstrained = findTransportRxns(male_unconstrained, '[luSI]', '[luLI]', false);
luSI_transRxnsT_female_unconstrained = findTransportRxns(female_unconstrained, '[luSI]', '[luLI]', false);
%1400 (female) 1402 (male) rxns have 0 to 1000 bound

luC_transRxnsT_male_unconstrained = findTransportRxns(male_unconstrained, '[luC]', '[luLI]', false);
luC_transRxnsT_female_unconstrained = findTransportRxns(female_unconstrained, '[luC]', '[luLI]', false);
%376 (female) 375 (male rxns have 0 to 1000 bound

fe_transRxnsT_male_unconstrained = findTransportRxns(male_unconstrained, '[fe]', '[luLI]', false);
fe_transRxnsT_female_unconstrained = findTransportRxns(female_unconstrained, '[fe]', '[luLI]', false);
%1398 (female) 1400 (male rxns have 0 to 1000 bound


load('LIECfemale.mat')
[femboundsDiff_LIEC, femdictionaryMap_LIEC, femuniquePathways_LIEC, femuniqueMets_LIEC] = compareStoichiometry(LIECfemale, subModel_female_Consts, 'Colon', 'Sub',0);



[maletransportSummary, maletransportMetsList, maleallRxns] = analyzeTransportsAndExchanges(colonModelMale_withDietFe);
[femaletransportSummary, femaletransportMetsList, femaleallRxns] = analyzeTransportsAndExchanges(colonModelFemale_withDietFe);

[subfemaletransportSummary, subfemaletransportMetsList, subfemaleallRxns] = analyzeTransportsAndExchanges(subModel_female_Consts);
[submaletransportSummary, submaletransportMetsList, submaleallRxns] = analyzeTransportsAndExchanges(subModel_male_Consts);


[femboundsDiff_withDietFe, femdictionaryMap_withDietFe, femuniquePathways_withDietFe, femuniqueMets_withDietFe] = compareStoichiometry(colonModelMale_withDietFe, subModel_female_Consts, 'Colon', 'Sub',0);
[maleboundsDiff_withDietFe, maledictionaryMap_withDietFe, maleuniquePathways_withDietFe, maleuniqueMets_withDietFe] = compareStoichiometry(colonModelMale_withDietFe, subModel_male_Consts, 'Colon', 'Sub',0);

[femboundsDiff, femdictionaryMap, femuniquePathways, femuniqueMets] = compareStoichiometry(ColonModelFemale, subModel_female_Consts, 'Colon', 'Sub',0);
[maleboundsDiff, maledictionaryMap, maleuniquePathways, maleuniqueMets] = compareStoichiometry(ColonModelMale, subModel_male_Consts, 'Colon', 'Sub',0);

% 1. Define the name of the Excel file you want to create
filename = fullfile(resultsPath, 'BoundaryDifferences_by_Category_male.xlsx');

% 2. Extract all unique categories from your table
maleuniqueCats = unique(maleboundsDiff_withDietFe.Category);

fprintf('\nPreparing to save %d categories to Excel...\n', length(maleuniqueCats));

% 3. Loop through each category and save it to its own sheet
for i = 1:length(maleuniqueCats)
    currentCat = maleuniqueCats{i};
    
    % Filter the master table for only the current category
    catTable = maleboundsDiff_withDietFe(strcmp(maleboundsDiff_withDietFe.Category, currentCat), :);
    
    % Excel Safety Check: Sheet names cannot be longer than 31 characters
    sheetName = currentCat;
    if length(sheetName) > 31
        sheetName = sheetName(1:31);
    end
    
    % Write the filtered table to a specific sheet in the Excel file
    %writetable(catTable, filename, 'Sheet', sheetName);
    
    fprintf(' -> Saved %3d rxns to sheet: %s\n', height(catTable), sheetName);
end

fprintf('Successfully saved all categories to %s\n', filename);


% 1. Define the name of the Excel file you want to create
filename = fullfile(resultsPath, 'BoundaryDifferences_by_Category_Female.xlsx');

% 2. Extract all unique categories from your table
femuniqueCats = unique(femboundsDiff_withDietFe.Category);

fprintf('\nPreparing to save %d categories to Excel...\n', length(femuniqueCats));

% 3. Loop through each category and save it to its own sheet
for i = 1:length(femuniqueCats)
    currentCat = femuniqueCats{i};
    
    % Filter the master table for only the current category
    catTable = femboundsDiff_withDietFe(strcmp(femboundsDiff_withDietFe.Category, currentCat), :);
    
    % Excel Safety Check: Sheet names cannot be longer than 31 characters
    sheetName = currentCat;
    if length(sheetName) > 31
        sheetName = sheetName(1:31);
    end
    
    % Write the filtered table to a specific sheet in the Excel file
    %writetable(catTable, filename, 'Sheet', sheetName);
    
    fprintf(' -> Saved %3d rxns to sheet: %s\n', height(catTable), sheetName);
end

fprintf('Successfully saved all categories to %s\n', filename);

