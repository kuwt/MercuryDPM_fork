function CDF = convertCDFtoPDF(PDF)
    CDF = zeros(length(PDF),2);
    for i = 2:length(PDF)
        PDF(i,2) = min(1, PDF(i,2) + PDF(i-1,2));
    end
    CDF = PDF;
    CDF = validateCDF(CDF);