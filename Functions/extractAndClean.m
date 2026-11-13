function T = extractAndClean(model, rxnIdx, prefix, cleanFunc)
    rxns = model.rxns(rxnIdx);
    lbs = model.lb(rxnIdx);
    ubs = model.ub(rxnIdx);
    
    baseMets = cleanFunc(rxns);
    T = table(baseMets, lbs, ubs, 'VariableNames', {'BaseMetabolite', [prefix '_lb'], [prefix '_ub']});
    
    % Ensure unique BaseMetabolites to prevent outerjoin explosion
    [~, uniqueIdx] = unique(T.BaseMetabolite, 'stable');
    T = T(uniqueIdx, :);
end
