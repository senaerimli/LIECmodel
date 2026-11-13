function [moreThanOne, equalsToZero, lessThanOne] = calculateShadowPrice(model)
    % Run optimization
    sol = optimizeCbModel(model);
    
    % Create Logical Masks (1s and 0s)
    idxMore = sol.y > 1 | sol.y < -1;
    idxZero = sol.y == 0;
    
    % Everything else falls into the "less than one" category
    idxLess = ~(idxMore | idxZero); 
    
    % Extract data directly using the masks
    moreThanOne.mets = model.mets(idxMore);
    moreThanOne.y    = sol.y(idxMore);
    
    equalsToZero.mets = model.mets(idxZero);
    equalsToZero.y    = sol.y(idxZero);
    
    lessThanOne.mets = model.mets(idxLess);
    lessThanOne.y    = sol.y(idxLess);
end