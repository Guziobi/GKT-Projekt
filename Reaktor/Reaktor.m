% Reaktor KAA146 Grundläggande kemiteknik, Projektarbete grupp 4

clc, clear

% isobutan = A
% isobuten = B
% vätgas   = C
% vatten   = D

% Cp beräkning
% H2O
clc
clear
T=1100;
A_H2O = 30.092;
B_H2O = 6.832514;
C_H2O = 6.793435;
D_H2O = -2.53448;87
Cp_H2O = (A_H2O + (B_H2O.*(T.*10^-3)) + (C_H2O.*((T.*10^-3).^2)) + (D_H2O.*((T.*10^-3).^3)))%[kJ/mol]
% H2
A_H2 = 33.066178;
B_H2 = -11.363417;
C_H2 = 11.432816;
D_H2 = 2.772874;
Cp_H2 = (A_H2 + (B_H2.*(T.*10^-3)) + (C_H2.*((T.*10^-3).^2)) + (D_H2.*((T.*10^-3).^3)))

%
%Modell för reaktorn

% isobutan = A
% isobuten = B
% vätgas   = C
% vatten   = D

clc, clear

% Data
R = 8.31447;
Cp_D = 100;
l = [0 10];     % Längd på reaktor1

