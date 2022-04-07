% Reaktor KAA146 Grundläggande kemiteknik, Projektarbete grupp 4

clc, clear

% isobutan = A
% isobuten = B
% vätgas   = C
% vatten   = D

% Data
R = 8.31447;
rho_cat = 1120; % [kg/m3]
V = 10; %[m3]
Tin_reak = 750;
Ea = 141e3; %[J/mol]
k = 0.0596; % mol/kg cat.*s*bar vid 550 C
K1 = 22.9; % bar^-1
K2 = 7.56; % bar^-1

CPcoeff_H2O = [72.43 1.039*10^-2 -1.497*10^-6 0 ];
CPcoeff_H2 = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];
CPcoeff_ISOBUTAN = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];
CPcoeff_ISOBUTEN = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];

% Arrenius ekv.
A = k/exp(-Ea/(823*R));

% REAKTOR 1

% Cp beräkning
T = Tin_reak; %K, equlibrium-conversion mot T börjar avta vid denna ~temp

CpH2O = Cp_calc(T,CPcoeff_H2O);
CpH2 = Cp_calc(T,CPcoeff_H2);
CpButan = Cp_calc(T,CPcoeff_ISOBUTAN);
CpButen = Cp_calc(T,CPcoeff_ISOBUTEN);

Cp = [CpButan CpButen CpH2 CpH2O];

dH0_Butan   = -17.9*10^3;
dH0_Buten   = -134.2*10^3;
dH0_H2      = 0;

dH0         = dH0_Buten+dH0_H2-dH0_Butan;

            %Butan Buten H2 H2O Temp
U0          = [128 5 0 1091 750];

Vstart = 0; %Volym m3
Vfinal = V; 
Vspan = [Vstart Vfinal];
[V,U] = ode15s(@PFR_ode,[Vstart Vfinal],U0,[],Cp,dH0,A,Ea);

% REAKTOR 2

% Cp beräkning
T = Tin_reak; %K, equlibrium-conversion mot T börjar avta vid denna ~temp


plot(V, U(:,1))
