
% --- Helper Function for Clean Excel Exporting ---
function exportShadowPrices(priceStruct, fieldName, filePath, sheetName)
    if isfield(priceStruct, fieldName) && isfield(priceStruct.(fieldName), 'mets') && ~isempty(priceStruct.(fieldName).mets)
        T = struct2table(priceStruct.(fieldName));
        writetable(T, filePath, 'Sheet', sheetName);
    end
end