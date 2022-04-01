    % Cp beräkning
    function Cp = Cp_calc(T)
    % H2O
    A_H2O = 30.092;
    B_H2O = 6.832514;
    C_H2O = 6.793435;
    D_H2O = -2.53448;
    Cp_H2O = (A_H2O + (B_H2O.*(T.*10^-3)) + (C_H2O.*((T.*10^-3).^2)) + (D_H2O.*((T.*10^-3).^3)));%[kJ/mol]
    
    % H2
    A_H2 = 33.066178;
    B_H2 = -11.363417;
    C_H2 = 11.432816;
    D_H2 = 2.772874;
    Cp_H2 = (A_H2 + (B_H2.*(T.*10^-3)) + (C_H2.*((T.*10^-3).^2)) + (D_H2.*((T.*10^-3).^3)));%[kJ/mol]
    
    % Isobutan
    A_butan = -1.39;
    B_butan = 0.3847;
    C_butan = -1.846*10^-4;
    D_butan = 2.895*10^-8;
    Cp_butan = A_butan+(B_butan.*T)+(C_butan.*(T^2))+(D_butan.*(T^3));%[kJ/mol]
    
    % Isobuten
    A_buten = 16.05;
    B_buten = 0.2804;
    C_buten = -1.091*10^-4;
    D_buten = 9.098*10^-9;
    Cp_buten = A_buten+(B_buten.*T)+(C_buten.*(T^2))+(D_buten.*(T^3));%[kJ/mol]
    
    Cp = [Cp_butan Cp_buten Cp_H2 Cp_H2O]; %[A(butan) B(buten) C(H2) D(H2O)]
    end