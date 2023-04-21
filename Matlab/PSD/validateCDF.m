function [validatedCDF] = validateCDF(CDF)

% Ensure range of probability to be [0,1]
for i = 1:length(CDF)
    CDF(i,2) = CDF(i,2)/CDF(end,2);
end

% insert zero at the beginning
if CDF(1,2) ~= 0
    rowToInsert = [CDF(1,2),0];
    validatedCDF = [rowToInsert;CDF];
end
% insert one at beginning
if CDF(end,2) ~= 0
    rowToInsert = [CDF(end,2), 1];
    validatedCDF = [CDF;rowToInsert];
end
validatedCDF = CDF;
end