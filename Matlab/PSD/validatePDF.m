function [validatedPDF] = validatePDF(PDF)

% insert zero at the beginning
if PDF(1,2) ~= 0
%     rowtToInsert = [max(1 *PDF(1,1), 2 * PDF(1,1)-PDF(2,1)), 0];
    rowToInsert = [PDF(1,1),0];
    validatedPDF = [rowToInsert;PDF];
end
validatedPDF = PDF;
end