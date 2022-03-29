% Reaktor KAA146 Grundläggande kemiteknik, Projektarbete grupp 4

clc, clear

% isobutan = A
% isobuten = B
% vätgas   = C
% vatten   = D

% Data
R = 8.31447;
rho_cat = 1120; % [kg/m3]
V = 1; %[m3]
T_reaktorer = [298 832 298 832]; %K [T1 T2 T3 T4]
Ea = 141e3; %[J/mol]
k = 0.0596; % mol/kg cat.*s*bar vid 550 C
K1 = 22.9; % bar^-1
K2 = 7.56; % bar^-1

% REAKTOR 1

% Cp beräkning
T = T_reaktorer(2); %K, equlibrium-conversion mot T börjar avta vid denna ~temp
% H2O
A_H2O = 30.092;
B_H2O = 6.832514;
C_H2O = 6.793435;
D_H2O = -2.53448;
E_H2O = 0.082139;
Cp_H2O = A_H2O + (B_H2O.*(T.*10^-3)) + (C_H2O.*((T.*10^-3).^2)) + (D_H2O.*((T.*10^-3).^3))+(E_H2O./((T.*10^-3)^2));%[kJ/mol]

% H2
A_H2 = 33.066178;
B_H2 = -11.363417;
C_H2 = 11.432816;
D_H2 = 2.772874;
E_H2 = -0.158558;
Cp_H2 = A_H2 + (B_H2.*(T.*10^-3)) + (C_H2.*((T.*10^-3).^2)) + (D_H2.*((T.*10^-3).^3))+(E_H2./((T.*10^-3)^2));%[kJ/mol]

% REAKTOR 2

% Cp beräkning
T = T_reaktorer(4); %K, equlibrium-conversion mot T börjar avta vid denna ~temp
%H2O
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

