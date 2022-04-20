% Reaktor KAA146 Grundläggande kemiteknik, Projektarbete grupp 4

clc, clear, close all

% isobutan = A
% isobuten = B
% vätgas   = C
% vatten   = D

% Data
R = 8.31447;
P = 1;              % [bar]
T_reaktor1 = 1000;   % [K]
T_reaktor2 = 1000;   % [K]
rho_cat = 1120;     % [kg m3^-1]
W = 7850;           % [kg]
Ea = 141e3;         % [J mol^-1]
k = 0.0596;         % [mol kg cat.^-1 s^-1 bar^-1 vid 550 C]
K1 = 22.9;          % [bar^-1]
K2 = 7.56;          % [bar^-1]

CPcoeffA = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];     % ISOBUTAN
CPcoeffB = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];     % ISOBUTEN
CPcoeffC = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];    % H2
CPcoeffD = [72.43 1.039*10^-2 -1.497*10^-6 0 ];         % H2O

Cp = [CPcoeffA; CPcoeffB; CPcoeffC; CPcoeffD];

dH0A   = -134.2e3;       % [kJ mol^-1]
dH0B   = -17.9e3;        % [kJ mol^-1]
dH0C   = 0;              % [kJ mol^-1]
dHr0   = dH0B+dH0C-dH0A; % [kJ mol^-1]

% REAKTOR 1
U01 = [168/3.6 6/3.6 0 1680/3.6 T_reaktor1]; %[mol s^-1]

Wstart1 = 0; %Massa cat. [kg]
Wfinal1 = W; 
Wspan1 = [Wstart1 Wfinal1];
[W1,U1] = ode15s(@PFR_ode,Wspan1,U01,[],dHr0,k,K1,K2,P,Cp);

X1 = (U1(1,1)-U1(:,1))/U1(1,1);
T1 = U1(:,5); %[K]

V1= W1./rho_cat; %[m3]

figure(1);
plot(W1,X1)
title('Reaktor 1'), xlabel('W [kg]'), ylabel('Omsättningsgrad, X')

figure (2);
plot(T1,X1);
title('Reaktor 1'), xlabel('T [K]'), ylabel('Omsättningsgrad, X')


% REAKTOR 2
U02 = [U1(1,1) U1(1,2) U1(1,3) U1(1,4) T_reaktor2];
Wstart2 = 0; %Massa cat. [kg]
Wfinal2 = 7850; 
Wspan2 = [Wstart2 Wfinal2];
[W2,U2] = ode15s(@PFR_ode,Wspan2,U02,[],dHr0,k,K1,K2,P,Cp);

X2 = (U2(1,1)-U2(:,1))/U2(1,1);
T2 = U2(:,5); %[K]

V2= W2./rho_cat; %[m3]

figure(3);
plot(W2,X2)
title('Reaktor 2'),xlabel('W [kg]'), ylabel('Omsättningsgrad, X')

figure (4);
plot(T2,X2);
title('Reaktor 2'), xlabel('T [K]'), ylabel('Omsättningsgrad, X')

% BÅDA REAKTORERNA
W = [W1; [W2+W1(end)]]; % [kg]
X = [X1; [X2+X1(end)]]; 
T = [T1; T2];           % [K]
V = [V1; [V2+V1(end)]]; % [m3]

figure(5);
plot(W,X)
title('Reaktor 1 och 2'),xlabel('W [kg]'), ylabel('Omsättningsgrad, X')

figure(6);
plot(X,T)
title('Reaktor 1 och 2'),xlabel('Omsättningsgrad, X'),ylabel('T [K]')

figure(7);
plot(V,X)
title('Reaktor 1 och 2'),xlabel('V [m3]'), ylabel('Omsättningsgrad, X')

