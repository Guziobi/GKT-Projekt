function Cp = Cp_calc(T, CPcoeff)
    A = CPcoeff(1);
    B = CPcoeff(2);
    C = CPcoeff(3);
    D = CPcoeff(4);

    Cp = A + B*T + C*T.^2 + D*T.^3;

end

% CPcoeff_H2O = [72.43 1.039*10^-2 -1.497*10^-6 0 ];
% CPcoeff_H2 = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];
% CPcoeff_ISOBUTAN = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];
% CPcoeff_ISOBUTEN = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];