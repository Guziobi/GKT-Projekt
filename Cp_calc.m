function Cp = Cp_calc(T, CPcoeff)
    A = CPcoeff(1);
    B = CPcoeff(2);
    C = CPcoeff(3);
    D = CPcoeff(4);

    Cp = A + B*T + C*T^2 + D*T^3;

end

% CPcoeff_H2O = [30.092 6.832514 6.793435 -2.53448];
% CPcoeff_H2 = [33.066178 -11.363417 11.432816 2.772874];
% CPcoeff_ISOBUTAN = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];
% CPcoeff_ISOBUTEN = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];