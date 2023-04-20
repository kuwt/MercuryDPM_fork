function PNDF = convertPVDFtoPNDF(PVDF)

PNDF = PVDF;
k_v = 4/3 * pi;
for i = 2:length(PVDF)
    PNDF(i,2) = PVDF(i,2)/(k_v*0.25*((PVDF(i,1)^2+PVDF(i-1,1)^2)) * (PVDF(i,1) + PVDF(i,1)));
end
% PNDF(:,2) = PNDF(:,2)/sum(PNDF(:,2));
validatePDF(PNDF)
end