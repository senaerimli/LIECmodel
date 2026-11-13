function [mismatchedBoundsTable] = findBoundDiff(model1, model2)
% compareModelBoundsByS Compares the bounds of two models based on their S-matrix.
% Matches reactions using stoichiometry for shared metabolites, accounting
% for both forward and reverse reaction definitions. Filters output to 
% only show mismatched constraints.
%
% INPUTS:
%   model1 - First COBRA model structure
%   model2 - Second COBRA model structure
%
% OUTPUTS:
%   mismatchedBoundsTable - A MATLAB table containing the reactions 
%                           with matching stoichiometry but differing bounds.
%
% EXAMPLE:
%   mismatchedRxns = findBoundDiff(ColonModelFemale, modelWBM)
% Sena Erimli, 2026

% 1. Align the Metabolites (Rows)
   
[commonMets, idx1, idx2] = intersect(model1.mets, model2.mets);
    
if isempty(commonMets)
    error('No common metabolites found between the two models. Cannot compare S matrices.');
end

% Standardize "Infinite" Bounds
% Replaces magnitudes of 1,000 with 1,000,000 
model1.lb(model1.lb == -1000) = -1000000;
model1.ub(model1.ub ==  1000) =  1000000;

model2.lb(model2.lb == -1000) = -1000000;
model2.ub(model2.ub ==  1000) =  1000000;    

% 2. Extract and align S matrices 
% We transpose (') them immediately because ismember(..., 'rows') compares rows, 
% but our reactions are stored as columns.
% Note: full() is used because ismember requires full matrices for row-wise comparison.
    
S1_aligned = full(model1.S(idx1, :))';
S2_aligned = full(model2.S(idx2, :))';
    
% 3. Find Stoichiometric Matches (Forward and Reverse)
[isMatch_Fwd, loc_Fwd] = ismember(S1_aligned, S2_aligned, 'rows');
[isMatch_Rev, loc_Rev] = ismember(-S1_aligned, S2_aligned, 'rows');
    
% 4. Compare Bounds
mismatched_data = {};

for i = 1:length(model1.rxns) 
    rxn1 = model1.rxns{i};
    lb1 = model1.lb(i);
    ub1 = model1.ub(i);
        
    if isMatch_Fwd(i)
        idx2_match = loc_Fwd(i);
        rxn2 = model2.rxns{idx2_match};
        lb2 = model2.lb(idx2_match);
        ub2 = model2.ub(idx2_match);
             
        % If bounds are NOT identical, record them
        if (lb1 ~= lb2) || (ub1 ~= ub2)
            mismatched_data(end+1, :) = {rxn1, rxn2, 'Exact Match', lb1, lb2, ub1, ub2}; 
        end
            
    elseif isMatch_Rev(i)
        idx2_match = loc_Rev(i);
        rxn2 = model2.rxns{idx2_match};
        lb2 = model2.lb(idx2_match);
        ub2 = model2.ub(idx2_match);
            
        % For reverse matches, M1's lb maps to M2's -ub, and M1's ub maps to M2's -lb

        if (lb1 ~= -ub2) || (ub1 ~= -lb2)
            mismatched_data(end+1, :) = {rxn1, rxn2, 'Reverse Match', lb1, lb2, ub1, ub2}; 
        end
    end
end
    
% 5. Format Output as a Table
if isempty(mismatched_data)
    fprintf('All %d stoichiometrically matched reactions have identical bounds.\n', sum(isMatch_Fwd | isMatch_Rev));
    mismatchedBoundsTable = table(); % Return an empty table
else
    mismatchedBoundsTable = cell2table(mismatched_data, ...
        'VariableNames', {'Model1_Rxn', 'Model2_Rxn', 'MatchType', ...
        'M1_lb', 'M2_lb', 'M1_ub', 'M2_ub'});
    fprintf('Found %d reactions with matching stoichiometry but mismatched bounds.\n', height(mismatchedBoundsTable));
end
end