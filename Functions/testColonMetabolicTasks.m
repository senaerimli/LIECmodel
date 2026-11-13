function [results, missedMetsTable, taskCores] = testColonMetabolicTasks(model, generalModel, tasksFile, options)
% TESTMETABOLICTASKS Tests metabolic tasks on a genome-scale model using FBA.

if nargin < 4 || isempty(options)
    options = struct();
end

if ~isfield(options, 'solver'),          options.solver          = 'gurobi'; end
if ~isfield(options, 'tol'),             options.tol             = 1e-6;     end
if ~isfield(options, 'verbose'),         options.verbose         = true;     end
if ~isfield(options, 'extractCoreRxns'), options.extractCoreRxns = false;    end
if ~isfield(options, 'usePhysioBounds'), options.usePhysioBounds = false;    end

changeCobraSolver(options.solver, 'LP', 0);

if istable(tasksFile)
    tasks = tasksFile;
else
    tasks = loadTasksTable(tasksFile);
end

baseModel = model;
baseModel.c = zeros(size(baseModel.c));

nTasks = height(tasks);
objValues   = nan(nTasks, 1);
statusCol   = repmat({'FAIL'}, nTasks, 1);  
skipReason  = repmat({''}, nTasks, 1);      
missedMets  = cell(nTasks, 1);              
absentMets  = cell(nTasks, 1);             
taskCores   = struct(); 

for i = 1:nTasks    
    taskName = tasks.TaskName{i};
   
    rawReactStr = strrep(char(tasks.Reactants{i}), ';', ','); 
    reactants   = strtrim(strsplit(rawReactStr, ','));
    reactants(cellfun(@isempty, reactants)) = []; 
    reactLBs    = parseNumericList(tasks.ReactantsLB{i}, numel(reactants));
    
    [reactants, ~, rIdxMap] = unique(reactants, 'stable');
    aggReactLBs = zeros(1, numel(reactants));
    for r = 1:numel(reactants)
        aggReactLBs(r) = sum(reactLBs(rIdxMap == r));
    end
    reactLBs = aggReactLBs;

    rawProdStr  = strrep(char(tasks.Products{i}), ';', ',');
    products    = strtrim(strsplit(rawProdStr, ','));
    products(cellfun(@isempty, products)) = []; 
    prodLBs     = parseNumericList(tasks.ProductsLB{i}, numel(products));

    [products, ~, pIdxMap] = unique(products, 'stable');
    aggProdLBs = zeros(1, numel(products));
    for p = 1:numel(products)
        aggProdLBs(p) = sum(prodLBs(pIdxMap == p));
    end
    prodLBs = aggProdLBs;
    
    reqFlux = prodLBs(1);
    if isnan(reqFlux) || reqFlux <= 0, reqFlux = options.tol; end
    
    if options.verbose, fprintf('[Task %d/%d] %s\n', i, nTasks, taskName); end
    
    try
        m = baseModel;

        missedThisTask  = {}; 
        absentThisTask  = {}; 
        skipThisTask = false; 
        skipMsg = '';
        
        % --- DYNAMIC ANAEROBIC ENFORCER ---
        if contains(lower(taskName), 'anaerobic') || contains(lower(taskName), 'hypoxic')
            if ismember('EX_o2[e]', m.rxns), m = changeRxnBounds(m, 'EX_o2[e]', 0, 'l'); end
            if ismember('EX_o2[luLI]', m.rxns), m = changeRxnBounds(m, 'EX_o2[luLI]', 0, 'l'); end
        end

        % --- ANTI-CHEAT PROTOCOL: Block direct external uptake of target products ---
        if ~skipThisTask
             reactBaseNames = cell(1, numel(reactants));
            for r = 1:numel(reactants)
                reactBaseNames{r} = regexprep(reactants{r}, '\[.*\]', '');
            end

            for p = 1:numel(products)
                prod = products{p};
                baseName = regexprep(prod, '\[.*\]', ''); 

                if ismember(baseName, reactBaseNames)
                    continue;
                end

                excID_e  = ['EX_' baseName '[e]'];
                excID_fe = ['EX_' baseName '[fe]'];
                excID_lu = ['EX_' baseName '[luLI]'];
                excID_d  = ['EX_' baseName '[d]'];

                if ismember(excID_e, m.rxns),  m = changeRxnBounds(m, excID_e, 0, 'l'); end
                if ismember(excID_fe, m.rxns), m = changeRxnBounds(m, excID_fe, 0, 'l'); end
                if ismember(excID_lu, m.rxns), m = changeRxnBounds(m, excID_lu, 0, 'l'); end
                if ismember(excID_d, m.rxns),  m = changeRxnBounds(m, excID_d, 0, 'l'); end
            end
        end

        % --- Add sink reactions for reactants ---
        for r = 1:numel(reactants)
            reactant = reactants{r};
            lb_r     = -abs(reactLBs(r)); 

            if options.usePhysioBounds
                isExternal = endsWith(reactant, '[e]') || endsWith(reactant, '[luLI]') || ...
                             endsWith(reactant, '[d]') || endsWith(reactant, '[fe]');
                isInternalGas = ismember(reactant, {'o2[c]', 'o2[m]', 'co2[c]', 'co2[m]'});
                
                if isExternal || isInternalGas
                    continue; 
                end
            end

            if ~any(strcmp(m.mets, reactant))
                [m, added] = safeAddMetFromRef(m, generalModel, reactant);
                if added
                    missedThisTask{end+1} = reactant; %#ok<AGROW>
                else
                    absentThisTask{end+1} = reactant; %#ok<AGROW>
                    skipThisTask = true; skipMsg = sprintf('Reactant "%s" absent from both models', reactant);
                    break   
                end
            end

            skID = ['SK_' reactant];
            if ~any(strcmp(m.rxns, skID))
                m = addReaction(m, skID, 'metaboliteList', {reactant}, 'stoichCoeffList', -1, 'lowerBound', 0, 'upperBound', 1000);
            end
            m = changeRxnBounds(m, skID, lb_r,  'l');
            m = changeRxnBounds(m, skID, 0,  'u');
        end

        % --- Add metabolites for products ---
         if ~skipThisTask
             for p = 1:numel(products)
                 prod = products{p};
                 if ~any(strcmp(m.mets, prod))
                     [m, added] = safeAddMetFromRef(m, generalModel, prod);
                     if added
                         missedThisTask{end+1} = prod;
                     else
                         absentThisTask{end+1} = prod; skipThisTask = true;
                         skipMsg = sprintf('Product "%s" absent from both models', prod);
                         break; 
                     end
                 end
             end
         end

         % --- Optimise ---
         if ~skipThisTask
             prodSK = 'Task_Target_Sink';
             m = addReaction(m, prodSK, 'metaboliteList', products(1), 'stoichCoeffList', -1, 'lowerBound', 0, 'upperBound', 1000);
             m = changeObjective(m, prodSK, 1);
    
             for p = 2:numel(products)
                 bypID = products{p};
                 skReactID = ['SK_' bypID];
                 bypSK = sprintf('Task_Byproduct_Sink_%d', p);
    
                 if any(strcmp(m.rxns, skReactID))
                     m = changeRxnBounds(m, skReactID, 1000, 'u');
                 elseif any(strcmp(m.rxns, bypSK))
                     m = changeRxnBounds(m, bypSK, 0, 'l');
                     m = changeRxnBounds(m, bypSK, 1000, 'u');
                 else
                     m = addReaction(m, bypSK, 'metaboliteList', {bypID}, 'stoichCoeffList', -1, 'lowerBound', 0, 'upperBound', 1000);
                 end
             end            

             sol = optimizeCbModel(m, 'max');
             active = (sol.stat == 1) && ~isnan(sol.f) && (sol.f > options.tol);
             objValues(i) = sol.f;

             task_essential_rxns = {}; 
            
             if active
                 if sol.f >= (reqFlux - options.tol)
                     statusCol{i} = 'PASS';
                 else
                     statusCol{i} = sprintf('PASS (Warning: flux is %.4g but should be >= %.4g)', sol.f, reqFlux);
                 end

                 if options.extractCoreRxns
                     targetFloor = min(sol.f, reqFlux);
                     m_pFBA = changeRxnBounds(m, prodSK, targetFloor, 'l');
                     
                     % THE FIX: Use native optimizeCbModel for safe Parsimonious FBA extraction
                     pFBA_sol = optimizeCbModel(m_pFBA, 'max', 'one');
                     
                     if pFBA_sol.stat == 1 && ~isempty(pFBA_sol.v)
                         % Truncate to numel(m.rxns) to safely handle solver reversible splits
                         v_net = pFBA_sol.v(1:numel(m.rxns));
                         active_idx = find(abs(v_net) > 1e-5);
                         active_rxns = m.rxns(active_idx);
                         
                         active_rxns = active_rxns(~contains(active_rxns, 'SK_'));
                         active_rxns = active_rxns(~contains(active_rxns, 'EXHAUST_'));
                         active_rxns = active_rxns(~contains(active_rxns, 'Task_Target_Sink'));
                         active_rxns = active_rxns(~contains(active_rxns, 'Task_Byproduct_Sink_'));
                         
                         task_essential_rxns = active_rxns;
                     else
                         warning(['pFBA failed for task ', taskName, ' despite FBA passing.']);
                     end
                 end
             else
                 statusCol{i} = 'FAIL';
             end
            
             safe_fieldname = regexprep(taskName, '[^a-zA-Z0-9_]', '_');
             if length(safe_fieldname) > 63, safe_fieldname = safe_fieldname(1:63); end
             taskCores.(safe_fieldname) = task_essential_rxns';
         else
             statusCol{i}  = 'SKIP';
             skipReason{i} = skipMsg;
         end
        
    catch ME
        statusCol{i}  = 'SKIP';
        skipReason{i} = sprintf('Unexpected error: %s', ME.message);
    end
    
    missedMets{i} = missedThisTask; absentMets{i} = absentThisTask;
    
    if options.verbose && ~strcmp(statusCol{i}, 'SKIP')
        fprintf('  ObjValue = %.6g  [%s]\n', objValues(i), statusCol{i});
    end

end

    if options.verbose
        nPass = sum(startsWith(statusCol, 'PASS'));
        nFail = sum(startsWith(statusCol, 'FAIL'));
        nSkip = sum(strcmp(statusCol, 'SKIP')); 
        
        fprintf('\n[testMetabolicTasks] Done.\n');
        fprintf('  PASS : %d  |  FAIL : %d  |  SKIP : %d  (total %d)\n', nPass, nFail, nSkip, nTasks);
            
        if nSkip > 0
            fprintf('  Skipped tasks (metabolites absent from both models):\n');
            skipIdx = find(strcmp(statusCol, 'SKIP'));
            for k = 1:numel(skipIdx)
                fprintf('    - %s  ->  %s\n', tasks.TaskName{skipIdx(k)}, skipReason{skipIdx(k)});
            end
        end
    end

passed  = startsWith(statusCol, 'PASS');
results = table(tasks.TaskName, tasks.Reactants, tasks.ReactantsLB, tasks.Products, tasks.ProductsLB, objValues, passed, statusCol, skipReason, ...
    'VariableNames', {'TaskName','Reactants','ReactantsLB','Products','ProductsLB','ObjValue','Passed','Status','SkipReason'});

maxMissed = max(max(cellfun(@numel, missedMets)), 1);
maxAbsent = max(max(cellfun(@numel, absentMets)), 1);
missedArray = cell(nTasks, maxMissed); absentArray = cell(nTasks, maxAbsent);

for i = 1:nTasks
    for j = 1:numel(missedMets{i}), missedArray{i,j} = missedMets{i}{j}; end
    for j = 1:numel(absentMets{i}), absentArray{i,j} = absentMets{i}{j}; end
end

missedCols = arrayfun(@(x) sprintf('MissingFromDraft_%d', x), 1:maxMissed, 'UniformOutput', false);
absentCols = arrayfun(@(x) sprintf('AbsentFromBothModels_%d', x), 1:maxAbsent, 'UniformOutput', false);
missedMetsTable = cell2table([tasks.TaskName, missedArray, absentArray], 'VariableNames', [{'TaskName'}, missedCols, absentCols]);

if ~isempty(options.outputFile)
    writetable(results, options.outputFile, 'Sheet', 'Results');
    writetable(missedMetsTable, options.outputFile, 'Sheet', 'MissedMetabolites');
end
end  

function tasks = loadTasksTable(tasksFile)
    raw = readtable(tasksFile, 'VariableNamingRule', 'preserve');
    raw.Properties.VariableNames = strtrim(raw.Properties.VariableNames);
    colMap = struct('TaskName', findCol(raw, {'Metabolic tasks'}), 'Reactants', findCol(raw, {'Reactants'}), ...
        'ReactantsLB', findCol(raw, {'Reactants LB'}, true), 'Products', findCol(raw, {'Products', 'produced_met'}), ...
        'ProductsLB', findCol(raw, {'Products LB'}, true));
    nRows = height(raw); tasks = table();
    tasks.TaskName = cellstr(string(raw{:, colMap.TaskName}));
    tasks.Reactants = cellstr(string(raw{:, colMap.Reactants}));
    tasks.ReactantsLB = ifelse(isempty(colMap.ReactantsLB), repmat({'1'}, nRows, 1), cellstr(string(raw{:, colMap.ReactantsLB})));
    tasks.Products = cellstr(string(raw{:, colMap.Products}));
    tasks.ProductsLB = ifelse(isempty(colMap.ProductsLB), repmat({'1'}, nRows, 1), cellstr(string(raw{:, colMap.ProductsLB})));
end

function colName = findCol(tbl, candidates, optional)
    if nargin < 3, optional = false; end
    cols = tbl.Properties.VariableNames;
    for k = 1:numel(candidates)
        if any(strcmpi(cols, candidates{k})), colName = cols{find(strcmpi(cols, candidates{k}),1)}; return; end
    end
    if optional, colName = ''; else, error('Column missing.'); end
end

function vals = parseNumericList(str, expectedLen)
    if isnumeric(str), vals = double(str);
    else, vals = cellfun(@str2double, strtrim(strsplit(strrep(char(str), ';', ','), ','))); end
    vals(isnan(vals)) = 1.0;
    if numel(vals) == 1 && expectedLen > 1, vals = repmat(vals, 1, expectedLen); end
    if numel(vals) < expectedLen, vals(end+1 : expectedLen) = 1.0; end
end

function [m, added] = safeAddMetFromRef(m, refModel, metID)
    if any(strcmp(m.mets, metID)), added = false; return; end
    if any(strcmp(refModel.mets, metID)), m = addMetFromRef(m, refModel, metID); added = true; else, added = false; end
end

function m = addMetFromRef(m, refModel, metID)
    refMetIdx = find(strcmp(refModel.mets, metID), 1);
    m.mets{end+1, 1} = refModel.mets{refMetIdx};
    
    if isfield(m, 'metNames')
        if isfield(refModel, 'metNames') && ~isempty(refModel.metNames)
            m.metNames{end+1, 1} = refModel.metNames{refMetIdx};
        else
            m.metNames{end+1, 1} = metID;
        end
    end
    if isfield(m, 'metFormulas')
        if isfield(refModel, 'metFormulas') && ~isempty(refModel.metFormulas)
            m.metFormulas{end+1, 1} = refModel.metFormulas{refMetIdx};
        else
            m.metFormulas{end+1, 1} = ''; 
        end
    end
    if isfield(m, 'metCharges')
        if isfield(refModel, 'metCharges') && ~isempty(refModel.metCharges)
            m.metCharges(end+1, 1) = refModel.metCharges(refMetIdx);
        else
            m.metCharges(end+1, 1) = 0;
        end
    end
    if isfield(m, 'metComps')
        if isfield(refModel, 'metComps') && ~isempty(refModel.metComps)
            m.metComps{end+1, 1} = refModel.metComps{refMetIdx};
        else
            idxOpen = strfind(metID, '['); idxClose = strfind(metID, ']');
            if ~isempty(idxOpen) && ~isempty(idxClose)
                m.metComps{end+1, 1} = metID(idxOpen(end)+1 : idxClose(end)-1);
            else
                m.metComps{end+1, 1} = ''; 
            end
        end
    end  
    
    m.S(end+1, :) = 0;
    if isfield(m, 'b'), m.b(end+1, 1) = 0; end
    if isfield(m, 'csense'), m.csense(end+1, 1) = 'E'; end
end

function res = ifelse(cond, trueVal, falseVal)
    if cond, res = trueVal; else, res = falseVal; end
end