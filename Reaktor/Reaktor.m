% Reaktor KAA146 Grundläggande kemiteknik, Projektarbete grupp 4

clc, clear

% isobutan = A
% isobuten = B
% vätgas   = C
% vatten   = D

% Data
R = 8.31447;
P = 1;              % [bar]
T = 800;            % [K]
rho_cat = 1120;     % [kg/m3]
W = 5000;           % [kg]
Ea = 141e3;         % [J/mol]
k = 0.0596;         % mol/kg cat.*s*bar vid 550 C
K1 = 22.9;          % bar^-1
K2 = 7.56;          % bar^-1

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
U01 = [140/3.6 6/3.6 0 1190/3.6 800]; %[mol s^-1]

Wstart1 = 0; %Volym m3
Wfinal1 = W; 
Wspan1 = [Wstart1 Wfinal1];
[W1,U1] = ode15s(@PFR_ode,Wspan1,U01,[],dHr0,k,K1,K2,P,Cp);

X1 = (U1(1,1)-U1(:,1))/U1(1,1);
T1 = U1(:,5);

V1= W1./rho_cat;

figure(1);
plot(V1,X1)
title('Reaktor 1'), xlabel('V [m3]'), ylabel('Omsättningsgrad, X')

figure (2);
plot(T1,X1);
title('Reaktor 1'), xlabel('T [K]'), ylabel('Omsättningsgrad, X')


% REAKTOR 2
U02 = [U1(1,1) U1(1,2) U1(1,3) U1(1,4) 800];
Wstart2 = 0; %Volym m3
Wfinal2 = 5000; 
Wspan2 = [Wstart2 Wfinal2];
[W2,U2] = ode15s(@PFR_ode,Wspan2,U02,[],dHr0,k,K1,K2,P,Cp);

X2 = (U2(1,1)-U2(:,1))/U2(1,1);
T2 = U2(:,5);

V2= W2./rho_cat;

figure(3);
plot(V2,X2)
title('Reaktor 2'),xlabel('V [m3]'), ylabel('Omsättningsgrad, X')

figure (4);
plot(T2,X2);
title('Reaktor 2'), xlabel('T [K]'), ylabel('Omsättningsgrad, X')

% KOSTNADSBERÄKNING REAKTOR 2
ReaktorV2 = Wfinal2./rho_cat; %[m3]

% BÅDA REAKTORERNA
V = [V1; [V2+V1(end)]]; 
X = [X1; [X2+X1(end)]]; 
T = [T1; T2]; 
figure(5);
plot(V,X)
title('Reaktor 1 och 2'),xlabel('V [m3]'), ylabel('Omsättningsgrad, X')

figure(6);
plot(X,T)
title('Reaktor 1 och 2'),xlabel('Omsättningsgrad, X'),ylabel('T [K]')