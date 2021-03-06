% Reaktor KAA146 Grundläggande kemiteknik, Projektarbete grupp 4

clc, clear, close all

% isobutan = A
% isobuten = B
% vätgas   = C
% vatten   = D

% Data
R = 8.31447;
P = 1;              % [bar]
T_reaktor1 = 950;   % [K]
T_reaktor2 = 950;   % [K]
rho_cat = 1120;     % [kg m3^-1]
Ea = 141e3;         % [J mol^-1]
k = 0.0596;         % [mol kg cat.^-1 s^-1 bar^-1 vid 550 C]
K1 = 22.9;          % [bar^-1]
K2 = 7.56;          % [bar^-1]

CPcoeffA = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];     % ISOBUTAN
CPcoeffB = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];     % ISOBUTEN
CPcoeffC = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];    % H2
CPcoeffD = [32.24 0.001924 1.055*10^-5 -3.596*10^-9];   % H2O

Cp = [CPcoeffA; CPcoeffB; CPcoeffC; CPcoeffD];

dH0A   = -134.2e3;       % [kJ mol^-1]
dH0B   = -17.9e3;        % [kJ mol^-1]
dH0C   = 0;              % [kJ mol^-1]
dHr0   = dH0B+dH0C-dH0A; % [kJ mol^-1]

% REAKTOR 1
U01 = [131/3.6 0/3.6 0 1310/3.6 T_reaktor1]; %[mol s^-1]

Wstart1 = 0; %Massa cat. [kg]
Wfinal1 = 16000; 
Wspan1 = [Wstart1 Wfinal1];
[W1,U1] = ode15s(@PBR_ode,Wspan1,U01,[],dHr0,k,K1,K2,P,Cp);

T1 = U1(:,5); %[K]
V1= W1./rho_cat; %[m3]

% REAKTOR 1 alternativ
U0_alt1 = [131/3.6 0/3.6 0 1310/3.6 T_reaktor1]; %[mol s^-1]

Wstart_alt1 = 0; %Massa cat. [kg]
Wfinal_alt1 = 35000; 
Wspan_alt1 = [Wstart_alt1 Wfinal_alt1];
[W_alt1,U_alt1] = ode15s(@PBR_ode,Wspan_alt1,U0_alt1,[],dHr0,k,K1,K2,P,Cp);

T_alt = U_alt1(:,5); %[K]
V_alt= W_alt1./rho_cat; %[m3]

X_alt1= (U_alt1(1,1)-U_alt1(:,1))/U_alt1(1,1);

% REAKTOR 2
U02 = [U1(end,1) U1(end,2) U1(end,3) U1(end,4) T_reaktor2];
Wstart2 = 0; %Massa cat. [kg]
Wfinal2 = 19000; 
Wspan2 = [Wstart2 Wfinal2];
[W2,U2] = ode15s(@PBR_ode,Wspan2,U02,[],dHr0,k,K1,K2,P,Cp);

T2 = U2(:,5); %[K]
V2 = W2./rho_cat; %[m3]

% BÅDA REAKTORERNA
U = [U1; U2([2:end],:)];
X = (U(1,1)-U(:,1))/U(1,1);
W = [W1; [W2([2:end],:)+W1(end)]]; % [kg]
T = [T1; T2([2:end],:)];           % [K]
V = [V1; [V2([2:end],:)+V1(end)]]; % [m3]

FA = U(:,1);
FB = U(:,2);
FC = U(:,3); 

% PLOTTAR:

% Plottar omsättningsgraden mot massan katalysator för de båda reaktorerna
figure(1);
plot(W,X)
hold on
plot(W_alt1,X_alt1,'r--')
title('Reaktor 1 och 2'),xlabel('Katalysatormassa, W [kg]'), ylabel('Omsättningsgrad, X')
hold on
plot([Wfinal1 Wfinal1],[0 1],'k--')     % Avgränsning mellan reaktor 1 och 2

% Plottar temperaturen mot omsättningsgraden för de båda reaktorerna
figure(2);
plot(X,T)
title('Reaktor 1 och 2'),xlabel('Omsättningsgrad, X'),ylabel('T [K]')

% Plottar molflöden för de olika specierna
figure(3);
plot(W,FA,W,FB)
hold on
plot([Wfinal1 Wfinal1],[0 40],'k--')   % Avgränsning mellan reaktor 1 och 2
title('Molflöde för specier vs katalysatormassa')
xlabel('Katalysatormassa, W [kg]'), ylabel('Molflöde [mol/s]')
legend('isobutan','isobuten och vätgas')

% Tryckkätlsväggens (reaktorväggens) tjocklek
% REKTOR 1
T_F1 = (9/5)*(T_reaktor1-273.15)+32; % Temperatur i Farenheit, om över 900F måste rostfritt stål användas, max 1500F
Vol1 = Wfinal1./rho_cat;
rad1 = (Vol1/(4*pi)).^(1/3); % Radie på reaktor 1 [m]
D1 = 2*rad1; % Diameter på reaktor 1 [m]
Smax = 74.5*10^6; % Maximalt tillåtna spänningen över 900F [N/m2]
P_konc = P*10^5*1.1; % Konstruktionstryck, 10% mer än arbetstryck [Pa]
wall1= ((P_konc*D1)./((2*Smax*1)-(1.2*P_konc))).*10^3; %[mm]
if D1<=2
    if wall1<7
        wall1=7;
    end
end
if D1>2 && D1<=2.5
    if wall1<9
        wall1=9;
    end
end
Mantel1 = 8*pi*(rad1^2);
V_wall1 = Mantel1*(wall1*10^-3); %[m3]
rho_wall = 8000; %[kg/m3]
mass_wall1 = V_wall1.*rho_wall; %[kg]
% Horisontell reaktor pris
a = 12800;
b = 73;
n = 0.85;
cost_reak1 = a+(b*(mass_wall1^n)); % Kostnad för reaktor 1 år 2010 (CEPCI = 532.9)
cost_reak1_cat = cost_reak1*1.5; % Kostnad för reaktor 1 + katalysator år 2010 (CEPCI = 532.9)
cost1_2020 = (cost_reak1_cat*(569/532.9))*9.99*4; % Kostnad för reaktor 1 + katalysator samt montering år 2020 i SEK
disp(['Kostnad för reaktor 1 år 2020 (SEK):   ',num2str(cost1_2020)]), disp(' ')

% REKTOR 1 alt.
T_Falt = (9/5)*(T_reaktor1-273.15)+32; % Temperatur i Farenheit, om över 900F måste rostfritt stål användas, max 1500F
Vol_alt = Wfinal_alt1./rho_cat;
rad_alt = (Vol_alt/(4*pi)).^(1/3); % Radie på reaktor 1 [m]
D_alt = 2*rad_alt; % Diameter på reaktor 1 [m]
Smax = 74.5*10^6; % Maximalt tillåtna spänningen över 900F [N/m2]
P_konc = P*10^5*1.1; % Konstruktionstryck, 10% mer än arbetstryck [Pa]
wall_alt= ((P_konc*D_alt)./((2*Smax*1)-(1.2*P_konc))).*10^3; %[mm]
if D_alt<=2
    if wall_alt<7
        wall_alt=7;
    end
end
if D_alt>2 && D_alt<=2.5
    if wall_alt<9
        wall_alt=9;
    end
end
if D_alt>2.5 && D_alt<=3
     if wall_alt<10
        wall_alt=10;
    end
end
Mantel_alt = 8*pi*(rad_alt^2);
V_wall_alt = Mantel_alt*(wall_alt*10^-3); %[m3]
rho_wall = 8000; %[kg/m3]
mass_wall_alt = V_wall_alt.*rho_wall; %[kg]
% Horisontell reaktor pris
a = 12800;
b = 73;
n = 0.85;
cost_reak_alt = a+(b*(mass_wall_alt^n)); % Kostnad för reaktor 1 år 2010 (CEPCI = 532.9)
cost_reak_alt_cat = cost_reak_alt*1.5; % Kostnad för reaktor 1 + katalysator år 2010 (CEPCI = 532.9)
cost_alt_2020 = (cost_reak_alt_cat*(569/532.9))*9.99*4; % Kostnad för reaktor 1 + katalysator samt montering år 2020 i SEK
disp(['Kostnad för alternativ reaktor 1 år 2020 (SEK):   ',num2str(cost_alt_2020)]), disp(' ')

% REAKTOR 2
T_F2 = (9/5)*(T_reaktor2-273.15)+32; % Temperatur i Farenheit, om över 900F måste rostfritt stål användas, max 1500F
Vol2 = Wfinal2./rho_cat;
rad2 = (Vol2/(4*pi)).^(1/3); % Radie på reaktor 1 [m]
D2 = 2*rad2; % Diameter på reaktor 1 [m]
Smax = 74.5*10^6; % Maximalt tillåtna spänningen över 900F [N/m2]
P_konc = P*10^5*1.1; % Konstruktionstryck, 10% mer än arbetstryck [Pa]
wall2= ((P_konc*D2)./((2*Smax*1)-(1.2*P_konc))).*10^3; %[mm]
if D2<=2
    if (wall2)<7
        wall2=7;
    end
end
if D2>2 && D2<=2.5
    if wall2<9
        wall2=9;
    end
end
Mantel2 = 8*pi*(rad2^2);
V_wall2 = Mantel2*(wall2.*10^-3); %[m3]
rho_wall = 8000; %[kg/m3]
mass_wall2 = V_wall2.*rho_wall; %[kg]
% Horisontell reaktor
a = 12800;
b = 73;
n = 0.85;
cost_reak2 = a+(b*(mass_wall2^n)); % Kostnad för reaktor 1 år 2010 (CEPCI = 532.9)
cost_reak2_cat = cost_reak2*1.5; % Kostnad för reaktor 1 + katalysator år 2010 (CEPCI = 532.9)
cost2_2020 = (cost_reak2_cat*(569/532.9))*9.99*4; % Kostnad för reaktor 1 + katalysator samt montering år 2020 i SEK
disp(['Kostnad för reaktor 2 år 2020 (SEK):   ',num2str(cost2_2020)]), disp(' ')
cost_allareakt = cost1_2020 + cost2_2020;
disp(['Kostnad för reaktor 1 och 2 år 2020 (SEK):   ',num2str(cost_allareakt)]), disp(' ')
%%
% Plottar 

function dUdV = PBR_ode(W, U, dHr0,k,K1,K2,P,Cp)

R = 8.31447;               % gaskon.

% A = isobutAn (C4H10)
% B = isobutEn (C4H8)
% C = Vätgas   (H2)
% D = vatten   (H2O)

% Hämtar molflöden och temperatur
FA = U(1); 
FB = U(2); 
FC = U(3); 
FD = U(4); 
T = U(5);
F_tot = FA + FB + FC + FD;

% Cp-funktioner beroende av T
CpA = @(T) Cp_calc(T,Cp(1,:));  %Cp för ämne A (isobutan) vid T
CpB = @(T) Cp_calc(T,Cp(2,:));  %Cp för ämne B (isobuten) vid T
CpC = @(T) Cp_calc(T,Cp(3,:));  %Cp för ämne C (H2) vid T
CpD = @(T) Cp_calc(T,Cp(4,:));  %Cp för ämne D (H2O) vid T

%Beräknar partialtryck
PA = P*(FA/F_tot);  %Partialtryck för ämne A (isobutan) för inflöde till reaktor
PB = P*(FB/F_tot);  %Partialtryck för ämne B (isobuten)) för inflöde till reaktor
PC = P*(FC/F_tot);  %Partialtryck för ämne C (H2)) för inflöde till reaktor
PD = P*(FD/F_tot);  %Partialtryck för ämne D (H2O) för inflöde till reaktor

% Ke som funktion av T
Ke = (2.1*10^7) * exp(-122*10^3/(8.314*T)); % [bar]

% Beräknar r
r = k*(PA - ((PB*PC)/Ke))/(1+K1*PB*((PC)^0.5)+(K2*PC)^0.5); % [mol kg cat.^-1 s^-1 bar^-0.5]

% Massbalanser för A -> B + C
dFA = -r;
dFB =  r;
dFC =  r;
dFD =  0;

% deltaH för A -> B + C
dHr = dHr0 + integral(@(T) CpA(T)-CpB(T)-CpC(T),298.15,T);  % Standardreaktionsentalpin beräknas vid 298.15 K

% Temperaturbalans
dTdV = r*(-dHr) / ((CpA(T)*FA + CpB(T)*FB + CpC(T)*FC + CpD(T)*FD));

dUdV = [dFA 
        dFB 
        dFC 
        dFD 
        dTdV];
end

function Cp = Cp_calc(T, CPcoeff)
    A = CPcoeff(1);
    B = CPcoeff(2);
    C = CPcoeff(3);
    D = CPcoeff(4);

    Cp = A + B*T + C*T.^2 + D*T.^3;

end