function vacuumedModel = vacuumTheModel(model)


% 1. Find all physiological exchanges (mirroring your prepareBaseModel)
all_boundaries = model.rxns(findExcRxns(model));

% 2. Apply your boundary setup:
model = changeRxnBounds(model, all_boundaries, 0, 'b');
        % 3. THIRD: Re-open the Wastes on 'm'
        % (Added Lactate and Alanine for anaerobic tasks)
        true_wastes = {'EX_co2[e]', 'EX_co2[fe]', 'EX_co2[luLI]', ...
                       'EX_hco3[e]', 'EX_hco3[fe]', 'EX_hco3[luLI]', ...
                       'EX_nh4[e]', 'EX_nh4[fe]', 'EX_nh4[luLI]', ...
                       'EX_pi[e]', 'EX_pi[fe]', 'EX_pi[luLI]', ...
                       'EX_urea[e]', 'EX_urea[fe]', 'EX_urea[luLI]', ...
                       'EX_so4[e]', 'EX_so4[fe]', 'EX_so4[luLI]', ...
                       'EX_lac_L[e]', 'EX_lac_L[fe]', 'EX_lac_L[luLI]', ...
                       'EX_lac_D[e]', 'EX_lac_D[fe]', 'EX_lac_D[luLI]', ...
                       'EX_ala_L[e]', 'EX_ala_L[fe]', 'EX_ala_L[luLI]'};
        valid_wastes = intersect(true_wastes, model.rxns);
        model = changeRxnBounds(model, valid_wastes, 1000, 'u');
        model = changeRxnBounds(model, valid_wastes, 0, 'l');
        
    fprintf('  [vacuumTheModel] Closed external exchange reactions.\n');
        
      % 3b. SPECIFIC BYPRODUCT EXHAUST (Strictly Excretion Only!)
        % ONLY physiological byproducts allowed. NO ENERGY COFACTORS.
        safe_vents = {
            'h2o2[c]', 'h2o2[m]', 'h2o2[x]', ... % Oxidases (MAO, DAO, Peroxisome)
            'no[c]', 'no[m]', ...                % Nitric Oxide gas
            'ahcys[c]', 'ahcys[m]', ...          % Methylation dead-end (Creatine, COMT)
            'lac_L[c]', 'lac_D[c]', 'ala_L[c]'   % Anaerobic fermentation wastes
        };
        
        for j = 1:length(safe_vents)
            met_name = safe_vents{j};
            if ismember(met_name, model.mets)
                clean_name = strrep(strrep(met_name, '[', '_'), ']', '');
                sinkID = ['EXHAUST_', clean_name];
                model = addReaction(model, sinkID, 'metaboliteList', {met_name}, 'stoichCoeffList', [-1], 'reversible', false);
                model = changeRxnBounds(model, sinkID, 1000, 'u');
                model = changeRxnBounds(model, sinkID, 0, 'l');
            end
        end
        
        %        
        % 4. FOURTH: Re-open the Solvents on 'm'
        solvents = {'EX_h2o[e]', 'EX_h2o[fe]', 'EX_h2o[luLI]', ...
                    'EX_h[e]', 'EX_h[fe]', 'EX_h[luLI]', ...
                    'EX_na1[e]', 'EX_na1[fe]', 'EX_na1[luLI]', ...
                    'EX_k[e]', 'EX_k[fe]', 'EX_k[luLI]', ...
                    'EX_cl[e]', 'EX_cl[fe]', 'EX_cl[luLI]'};
        valid_solvents = intersect(solvents, model.rxns);
        model = changeRxnBounds(model, valid_solvents, -1000, 'l');
        model = changeRxnBounds(model, valid_solvents, 1000, 'u');
        
        vacuumedModel = model;

end