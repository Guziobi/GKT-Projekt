function res = ekonomi(NTU, Cmin, Cmax, Ka, U, THin, TCin, tdrift, beta)

% syms NTU
% epsilon = (1 - exp(-NTU.*(1 - Cmin/Cmax)))/(1 - (Cmin/Cmax).*exp(-NTU.*(1 - Cmin/Cmax)));
% diff(epsilon);

dedNTU = (7840.*exp(-(49.*NTU)./209).*(exp(-(49.*NTU)./209) - 1)) ./(43681.*((160.*exp(-(49.*NTU)./209))./209 - 1).^2) - (49.*exp(-(49.*NTU)./209))./(209.*((160.*exp(-(49.*NTU)./209))./209 - 1));

bosse = Ka./(U.*(THin - TCin).*tdrift.*beta);

res = dedNTU - bosse;


end