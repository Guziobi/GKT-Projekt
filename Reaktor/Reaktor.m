% Reaktor KAA146 Grundläggande kemiteknik, Projektarbete grupp 4

clc, clear

% isobutan = A
% isobuten = B
% vätgas   = C
% vatten   = D

% Data
R = 8.31447;
<<<<<<< HEAD
P = 1;              % [bar]
rho_cat = 1120;     % [kg/m3]
W = 7000;              %[kg]
Ea = 141e3;         %[J/mol]
k = 0.0596;         % mol/kg cat.*s*bar vid 550 C
K1 = 22.9;          % bar^-1
K2 = 7.56;          % bar^-1

% REAKTOR 1
dH0_Butan   = -134.2*10^3;
dH0_Buten   = -17.9*10^3;
=======
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
>>>>>>> 96bd975e9251b4d2352214d5dba87c9bb60fd6f5
dH0_H2      = 0;

dH0 = dH0_Buten+dH0_H2-dH0_Butan;

<<<<<<< HEAD
U01 = [128/3.6 5/3.6 0 1091/3.6 750];

Wstart1 = 0; %Volym m3
Wfinal1 = W; 
Wspan1 = [Wstart1 Wfinal1];
[W1,U1] = ode15s(@PFR_ode,Wspan1,U01,[],dH0,k,K1,K2,P);

X1 = (U1(1,1)-U1(:,1))/U1(1,1);
T1 = U1(:,5);

figure(1);
plot(W1,X1)
title('Reaktor 1'), xlabel('W [kg]'), ylabel('Omsättningsgrad, X')

figure (2);
plot(T1,X1);
title('Reaktor 1'), xlabel('T [K]'), ylabel('Omsättningsgrad, X')
=======
            %Butan Buten H2 H2O Temp
U0          = [128 5 0 1091 750];
>>>>>>> 96bd975e9251b4d2352214d5dba87c9bb60fd6f5


% REAKTOR 2
U02 = [U1(1,1) U1(1,2) U1(1,3) U1(1,4) 750];
Wstart2 = 0; %Volym m3
Wfinal2 = 6000; 
Wspan2 = [Wstart2 Wfinal2];
[W2,U2] = ode15s(@PFR_ode,Wspan2,U02,[],dH0,k,K1,K2,P);

X2 = (U2(1,1)-U2(:,1))/U2(1,1);
T2 = U2(:,5);

figure(3);
plot(W2,X2)
title('Reaktor 2'),xlabel('W [kg]'), ylabel('Omsättningsgrad, X')

figure (4);
plot(T2,X2);
title('Reaktor 2'), xlabel('T [K]'), ylabel('Omsättningsgrad, X')

% BÅDA REAKTORERNA
W = [W1; [W2+W1(end)]]; 
X = [X1; [X2+X1(end)]]; 
T = [T1; T2]; 
figure(5);
plot(W,X)
title('Reaktor 1 och 2'),xlabel('W [kg]'), ylabel('Omsättningsgrad, X')

figure(6);
plot(X,T)
title('Reaktor 1 och 2'),xlabel('Omsättningsgrad, X'),ylabel('T [K]')

<<<<<<< HEAD
=======
% Cp beräkning
T = Tin_reak; %K, equlibrium-conversion mot T börjar avta vid denna ~temp


plot(V, U(:,1))
>>>>>>> 96bd975e9251b4d2352214d5dba87c9bb60fd6f5
