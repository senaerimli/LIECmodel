function updModel = updateColonSpecBoundaries(model)

%model = changeRxnBounds(model, 'ARGN', 0, 'b');  %it actually exists in
%healthy colonocyte
model = addReaction(model, 'SPMDSYN', 'reactionFormula', 'ametam[c] + ptrc[c] -> h[c] + 5mta[c] + spmd[c]');

%% Palmitoyl-CoA to the 3 R-group pools
model = addReaction(model, 'Map_pmt_to_Rtotal', ...
    'reactionFormula', 'pmtcoa[c] <=> Rtotalcoa[c]');

model = addReaction(model, 'Map_pmt_to_Rtotal2', ...
    'reactionFormula', 'pmtcoa[c] <=> Rtotal2coa[c]');

model = addReaction(model, 'Map_pmt_to_Rtotal3', ...
    'reactionFormula', 'pmtcoa[c] <=> Rtotal3coa[c]');

% 2. Ensure they are fully open
model = changeRxnBounds(model, {'Map_pmt_to_Rtotal', 'Map_pmt_to_Rtotal2', 'Map_pmt_to_Rtotal3'}, -1000, 'l');
model = changeRxnBounds(model, {'Map_pmt_to_Rtotal', 'Map_pmt_to_Rtotal2', 'Map_pmt_to_Rtotal3'}, 1000, 'u');

disp('Metabolite mapping bridges added!');

%% Phenylalanine to Tyrosine
model = changeRxnBounds(model, {'PHETHPTOX2','r0399'}, 0, 'b');

%% Bile acid synthesis - cholesterol metabolism
model = changeRxnBounds(model, {'RE1835C', 'RE1834C', 'RE1834M'}, 0, 'b');

%% Urea Cycle
model = changeRxnBounds(model, {'ARGSS', 'ARGSL'}, 0, 'b');
model = changeRxnBounds(model, 'RE2888C', 0, 'b'); %you need to delete this
model = changeRxnBounds(model, 'sink_arg_L[c]', 0, 'b');


%% Glutamine anaerobic oxidation to ATP
rxns_lb_0 = {'SUCCt4_2[luLI]'; 'SUCCt4_2[e]'; 'Transport_succ[luLI]'; ...
    'FUMtr[luLI]'; 'FUMtr[e]'; 'MAL_Lte[luLI]'; 'MAL_Lte[e]';...
    '3OHGLUTACOAT3t[e]'; 'HOMOVALte[e]'; '2HYDOGOAT3t[e]'; 'GLUTAROAT3t[e]';...
    'GLUTAROAT3t[luLI]'; 'r0879'; 'r2101'; 'r2102'; 'r2103'; 'r2104';...
    'r2105'; 'r2106'; 'r2107'; 'r2108'; 'r2109'; '2OXOADPt'; 'Transport_akg[luLI]';...
    'ABUTt2r[luLI]'; 'ABUTt4_2_r[luLI]'; 'ABUTt4_2_r[e]';'ABUTt2r';...
    'Transport_4abut[luLI]'};

model = changeRxnBounds(model, rxns_lb_0, 0, 'l');

rxns_ub_0 = {'FOLOAT1tc[luLI]'; 'CE1918te[e]'; 'FOLOAT1tc[e]'; 'KYNATEtr[luLI]';...
    'r2093'; 'KYNATEtr'; 'TMDOATthc'; 'SUCSALte'; '4OHBUTtce'};

model = changeRxnBounds(model, rxns_ub_0, 0, 'u');

disp('--- CLOSING ARTIFICIAL DEMAND SINKS ---');

% Find and block the GABA demand reaction (and any other potential amino acid/TCA demand artifacts)
% The 'b' sets both lower and upper bounds to 0.
demand_rxns_male = model.rxns(contains(model.rxns, 'DM_4abut') | contains(model.rxns, 'DM_sucsal'));
model = changeRxnBounds(model, demand_rxns_male, 0, 'b');
disp('Black holes closed!');

% 1. FIX THE PROPIONATE LEAK (Lock RE2649C to uptake only)
model = changeRxnBounds(model, 'RE2649C', 0, 'l');
model = changeRxnBounds(model, {'HMR_2088', 'HMR_2089'}, 0, 'l');

disp('--- ENFORCING H2O2 THERMODYNAMICS (MALE ONLY) ---');

% Find all H2O2, O2, and Superoxide metabolites in the male model
h2o2_mets = model.mets(startsWith(model.mets, 'h2o2['));
o2_mets = model.mets(startsWith(model.mets, 'o2['));
o2s_mets = model.mets(startsWith(model.mets, 'o2s['));

cheat_count = 0;

for i = 1:length(model.rxns)
    rxnID = model.rxns{i};
    
    % Skip boundary reactions safely
    if startsWith(rxnID, 'EX_') || startsWith(rxnID, 'DM_') || startsWith(rxnID, 'sink_')
        continue;
    end
    
    rxnMets = model.mets(model.S(:, i) ~= 0);
    h2o2_count = sum(ismember(rxnMets, h2o2_mets));
    
    % If exactly 1 H2O2 species is present
    if h2o2_count == 1
        has_o2_or_o2s = any(ismember(rxnMets, o2_mets)) || any(ismember(rxnMets, o2s_mets));
        
        % If it interacts with H2O2 WITHOUT O2 or Superoxide, it's illegal to run backward to synthesize it
        if ~has_o2_or_o2s
            for j = 1:length(h2o2_mets)
                idx = find(strcmp(model.mets, h2o2_mets{j}));
                if ~isempty(idx) && model.S(idx, i) ~= 0
                    coeff = model.S(idx, i);
                    
                    if coeff > 0 
                        % H2O2 is a product in the forward direction (->), so lock upper bound
                        model = changeRxnBounds(model, rxnID, 0, 'u');
                        cheat_count = cheat_count + 1;
                    elseif coeff < 0
                        % H2O2 is a reactant, meaning the reverse direction (<-) synthesizes it, so lock lower bound
                        model = changeRxnBounds(model, rxnID, 0, 'l');
                        cheat_count = cheat_count + 1;
                    end
                end
            end
        end
    end
end

disp(['Locked ', num2str(cheat_count), ' illegal peroxidases in model']);

disp('--- ENFORCING UNIVERSAL OXYGEN LAW (FIXED) ---');
o2_mets = model.mets(startsWith(model.mets, 'o2['));
cheat_count = 0;

for i = 1:length(model.rxns)
    rxnID = model.rxns{i};
    
    % 1. Skip boundary and sink reactions safely
    if startsWith(rxnID, 'EX_') || startsWith(rxnID, 'DM_') || startsWith(rxnID, 'sink_') || startsWith(rxnID, 'SK_')
        continue;
    end
    
    % 2. Skip Catalase (it is physiologically supposed to produce O2 from H2O2)
    if contains(rxnID, 'CAT')
        continue;
    end
    
    rxnMets = model.mets(model.S(:, i) ~= 0);
    
    if any(ismember(rxnMets, o2_mets))
        % 3. Check if it's a transport reaction (O2 is both consumed and produced)
        o2_coeffs = model.S(ismember(model.mets, o2_mets), i);
        if any(o2_coeffs > 0) && any(o2_coeffs < 0)
            continue; % It is transporting O2 between compartments, SKIP!
        end
        
        % 4. Only now, strictly enforce that standard internal enzymes cannot produce O2
        for j = 1:length(o2_mets)
            idx = find(strcmp(model.mets, o2_mets{j}));
            if ~isempty(idx) && model.S(idx, i) ~= 0
                coeff = model.S(idx, i);
                
                if coeff > 0 
                    model = changeRxnBounds(model, rxnID, 0, 'u');
                    cheat_count = cheat_count + 1;
                elseif coeff < 0
                    model = changeRxnBounds(model, rxnID, 0, 'l');
                    cheat_count = cheat_count + 1;
                end
            end
        end
    end
end
disp(['Locked ', num2str(cheat_count), ' illegal oxygen-producing enzymes.']);

disp('--- REPAIRING SEVERED PROPIONATE PATHWAY ---');
% 1. Get exact indices for the trapped metabolites
idx_ppcoa = find(strcmp(model.mets, 'ppcoa[m]'));
idx_mmcoa = find(strcmp(model.mets, 'mmcoa_S[m]'));

% 2. Find the irreversible reaction separating them
pcc_logical = (model.S(idx_ppcoa, :) ~= 0) & (model.S(idx_mmcoa, :) ~= 0);
broken_pcc_rxns = model.rxns(pcc_logical);

% 3. Open the lower bound to allow Propionyl-CoA to flow forward into the TCA cycle!
model = changeRxnBounds(model, broken_pcc_rxns, -1000, 'l');
disp(['Successfully reconnected the TCA cycle at: ', strjoin(broken_pcc_rxns, ', ')]);

disp('--- ENFORCING TISSUE SPECIFICITY (LIVER/KIDNEY LOCKDOWN) ---');

% 1. Lock Creatine Synthesis (Kidney/Liver specific)
% Block all reactions interacting with Guanidinoacetate (gudac)
idx_gudac = find(strcmp(model.mets, 'gudac[c]'));
if ~isempty(idx_gudac)
    rxns_gudac = model.rxns(model.S(idx_gudac, :) ~= 0);
    model = changeRxnBounds(model, rxns_gudac, 0, 'b');
    disp(['Locked ', num2str(length(rxns_gudac)), ' illegal Creatine synthesis reactions.']);
end

% 2. Lock Choline Oxidation to Glycine (Liver specific)
% Block all reactions interacting with Betaine Aldehyde (betald) and Betaine (glyb)
idx_betald_m = find(strcmp(model.mets, 'betald[m]'));
idx_betald_c = find(strcmp(model.mets, 'betald[c]'));
idx_glyb_c   = find(strcmp(model.mets, 'glyb[c]'));
idx_glyb_m   = find(strcmp(model.mets, 'glyb[m]'));

liver_choline_mets = [idx_betald_m; idx_betald_c; idx_glyb_c; idx_glyb_m];
liver_choline_mets(liver_choline_mets == 0) = []; % Remove empty indices if any

if ~isempty(liver_choline_mets)
    liver_choline_logical = false(length(model.rxns), 1);
    for m = 1:length(liver_choline_mets)
        liver_choline_logical = liver_choline_logical | (model.S(liver_choline_mets(m), :) ~= 0)';
    end
    rxns_liver_choline = model.rxns(liver_choline_logical);
    model = changeRxnBounds(model, rxns_liver_choline, 0, 'b');
    disp(['Locked ', num2str(length(rxns_liver_choline)), ' Liver choline oxidation reactions.']);
end

%   consistentIndices = fastcc(model, 1e-4, 1);
 %   allIndices = 1:length(model.rxns);
  %  inconsistentIndices = setdiff(allIndices, consistentIndices);
   % inconsistentRxnNames = model.rxns(inconsistentIndices);
    %disp(inconsistentRxnNames)

%model_consistent = convertToConsistent(model);

updModel = model;
end