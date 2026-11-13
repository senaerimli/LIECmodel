function [summaryTable, pathwayMets, allRxnsTable] = analyzeTransportsAndExchanges(model)
% ANALYZETRANSPORTSANDEXCHANGES Identifies boundary reactions and automatically
% aggregates intracellular compartments into 'COMBINED' macro-categories.
%
% OUTPUTS:
%   summaryTable : MATLAB Table showing counts of reactions and bounds.
%   pathwayMets  : Struct containing just the unique base metabolites.
%   allRxnsTable : Master table of every reaction (includes COMBINED tags).

    fprintf('\n======================================================\n');
    fprintf(' Running Global Transport & Exchange Analysis...\n');
    fprintf('======================================================\n');

    % Define which compartments are strictly "Intracellular"
    % You can add or remove tags here based on your model's conventions
    intraComps = {'c', 'm', 'r', 'x', 'g', 'l', 'n', 'i'}; 

    % 1. Normalize Extreme Bounds
    lb_fix = sum(model.lb == -1000000);
    ub_fix = sum(model.ub ==  1000000);
    model.lb(model.lb == -1000000) = -1000;
    model.ub(model.ub ==  1000000) =  1000;
    
    % 2. Extract Compartments
    metComps = cell(length(model.mets), 1);
    for i = 1:length(model.mets)
        metID = model.mets{i};
        tok = regexp(metID, '\[([a-zA-Z0-9_]+)\]$', 'tokens');
        if isempty(tok)
            tok = regexp(metID, '_([a-zA-Z0-9]+)$', 'tokens');
        end
        
        if ~isempty(tok)
            metComps{i} = tok{1}{1};
        else
            metComps{i} = 'unknown';
        end
    end

    % 3. Setup Trackers
    standard_bounds = [-1000, 0, 1000];
    tracker = containers.Map('KeyType', 'char', 'ValueType', 'any');

    masterCategories = {};
    masterRxnIDs = {};
    masterLBs = [];
    masterUBs = [];

    % 4. Loop through all reactions
    for i = 1:length(model.rxns)
        metIdx = find(model.S(:, i));
        rxnMets = model.mets(metIdx);
        rxnComps = unique(metComps(metIdx));
        
        isValidTarget = false;
        categoryKey = '';

        % -- Logic A: Exact Compartment Categories --
        if length(metIdx) == 1
            isValidTarget = true;
            categoryKey = sprintf('Exchange_from_%s', rxnComps{1});
        elseif length(rxnComps) == 2
            isValidTarget = true;
            sortedComps = sort(rxnComps);
            categoryKey = sprintf('%s_AND_%s', sortedComps{1}, sortedComps{2});
        end

        if isValidTarget
            baseMets = regexprep(rxnMets, '\[([a-zA-Z0-9_]+)\]$|_([a-zA-Z0-9]+)$', '');
            
            % Save Exact Data
            if ~isKey(tracker, categoryKey)
                tracker(categoryKey) = struct('Total', 0, 'Default', 0, 'Custom', 0, 'Mets', {{}});
            end
            stats = tracker(categoryKey);
            stats.Total = stats.Total + 1;
            if ismember(model.lb(i), standard_bounds) && ismember(model.ub(i), standard_bounds)
                stats.Default = stats.Default + 1;
            else
                stats.Custom = stats.Custom + 1;
            end
            stats.Mets = [stats.Mets; baseMets];
            tracker(categoryKey) = stats;
            
            masterCategories{end+1, 1} = categoryKey;
            masterRxnIDs{end+1, 1} = model.rxns{i};
            masterLBs(end+1, 1) = model.lb(i);
            masterUBs(end+1, 1) = model.ub(i);

            % -- Logic B: MACRO Grouping (Intracellular Aggregation) --
            groupedComps = rxnComps;
            hasIntra = false;
            for c_idx = 1:length(groupedComps)
                if ismember(groupedComps{c_idx}, intraComps)
                    groupedComps{c_idx} = 'intra';
                    hasIntra = true;
                end
            end
            groupedComps = unique(groupedComps); % Collapse duplicates like c and m
            
            if hasIntra
                macroKey = '';
                if length(metIdx) == 1
                    macroKey = sprintf('COMBINED_Exchange_from_%s', groupedComps{1});
                elseif length(rxnComps) == 2
                    if length(groupedComps) == 1
                        % If it collapsed to 1, it means both compartments were intracellular
                        macroKey = 'COMBINED_intra_AND_intra';
                    else
                        sortedMacro = sort(groupedComps);
                        macroKey = sprintf('COMBINED_%s_AND_%s', sortedMacro{1}, sortedMacro{2});
                    end
                end

                % Save Macro Data (Acts as a second entry for the same reaction)
                if ~isempty(macroKey)
                    if ~isKey(tracker, macroKey)
                        tracker(macroKey) = struct('Total', 0, 'Default', 0, 'Custom', 0, 'Mets', {{}});
                    end
                    mStats = tracker(macroKey);
                    mStats.Total = mStats.Total + 1;
                    if ismember(model.lb(i), standard_bounds) && ismember(model.ub(i), standard_bounds)
                        mStats.Default = mStats.Default + 1;
                    else
                        mStats.Custom = mStats.Custom + 1;
                    end
                    mStats.Mets = [mStats.Mets; baseMets];
                    tracker(macroKey) = mStats;

                    masterCategories{end+1, 1} = macroKey;
                    masterRxnIDs{end+1, 1} = model.rxns{i};
                    masterLBs(end+1, 1) = model.lb(i);
                    masterUBs(end+1, 1) = model.ub(i);
                end
            end
        end
    end

    % 5. Compile Outputs
    keys = tracker.keys;
    if isempty(keys)
        summaryTable = table(); pathwayMets = struct(); allRxnsTable = table();
        return;
    end
    
    Category = cell(length(keys), 1);
    Total_Rxns = zeros(length(keys), 1);
    Default_Bounds = zeros(length(keys), 1);
    Custom_Bounds = zeros(length(keys), 1);
    Unique_Metabolites_Count = zeros(length(keys), 1);
    pathwayMets = struct();

    for k = 1:length(keys)
        currentKey = keys{k};
        stats = tracker(currentKey);
        
        Category{k} = currentKey;
        Total_Rxns(k) = stats.Total;
        Default_Bounds(k) = stats.Default;
        Custom_Bounds(k) = stats.Custom;
        Unique_Metabolites_Count(k) = length(unique(stats.Mets));
        
        cleanField = regexprep(currentKey, '[^a-zA-Z0-9_]', ''); 
        pathwayMets.(cleanField) = unique(stats.Mets);
    end

    summaryTable = table(Category, Total_Rxns, Default_Bounds, Custom_Bounds, Unique_Metabolites_Count);
    
    % Sort table: Put COMBINED tags at the bottom, and sort by Total_Rxns within groups
    isCombined = startsWith(Category, 'COMBINED');
    [~, sortIdx] = sortrows([isCombined, -Total_Rxns]);
    summaryTable = summaryTable(sortIdx, :);

    fprintf('Generating formulas for the detailed master table...\n');
    masterFormulas = printRxnFormula(model, 'rxnAbbrList', masterRxnIDs, 'printFlag', false);
    
    allRxnsTable = table(masterCategories, masterRxnIDs, masterFormulas, masterLBs, masterUBs, ...
        'VariableNames', {'Pathway_Category', 'Reaction_ID', 'Formula', 'lb', 'ub'});

    fprintf('Successfully generated tables.\n');
    disp(summaryTable);
    fprintf('======================================================\n');
end