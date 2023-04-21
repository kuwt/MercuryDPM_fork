function PDF = convertCDFtoPDF(CDF)
    PDF = zeros(length(CDF),2);
    PDF(:,1) = CDF(:,1);
    probabilityOld = 0;
    for i = 1:length(CDF)
            probabilityCDF = CDF(i,2);
            PDF(i,2) = CDF(i,2) - probabilityOld;
            probabilityOld = probabilityCDF;
    end
    PDF = validatePDF(PDF);