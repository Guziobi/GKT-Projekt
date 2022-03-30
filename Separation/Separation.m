% Separation
clc, clear

% Data - separation 
q = 1;             % kokvarmt flöde
F = 100;           % kmol h-1
P = 760;           % mmHg (1 atm)
xf = 0.50;         % molbråk buten
xd = 0.95;         % destillatbråk 
xb = 0.05;         % bottenbråk
R = 10;              % återflödesförhållande
M1 = 56.1063;      % g mol-1
M2 = 58.1222;      % g mol-1

         % A      B       C
Ant1 =  [15.7564 2132.42 -33.15];  % buten
Ant2 =  [15.6782 2154.90 -34.42];  % butan 

W12 = 0.48584; 
W21 = 1.64637; 

rho_butenL = 1;          % kgm-3
rho_butanL = 1;          % kgm-3

liqsurften = 1;         % dyn cm-1

% Beräkning av flöden 
B = F*((xf-xd)/(xb-xd));
D = F-B;

L = R*D; 
V = L+D;

l = L + q*F; 
v = l-B;

% Återkokare 
[gamma1, gamma2] = wilson(xb,W12,W21);
Tstart = 273.15;
TB=fsolve(@(T)find_Tb(T,xb,gamma1,gamma2,Ant1,Ant2,P),Tstart);
Psat1 = antoine(TB, Ant1);
y0 = (gamma1*xb.*Psat1)/P;
x1 = (v/l)*y0 + (B/l)*xb;

% Avdrivardel 
x(1) = x1;
i = 0; 
y(1)=y0;

while x<xf
    i = i+1; 
    [gamma1, gamma2] = wilson(x(i),W12,W21);
    Tstart = TB;
    TB=fsolve(@(T)find_Tb(T,x(i),gamma1,gamma2,Ant1,Ant2,P),Tstart);
    Psat1 = antoine(TB, Ant1);                                                % Ångtryck
    y(i) = (gamma1*x(i).*Psat1)/P;
    x(i+1) = (v/l)*y(i) + (B/l)*xb;
end

m = i+1;

% Förstärkare
while y(i)<xd
    x(i+1)= (V/L)*y(i) + (1/L)*(B*xb-F*xf);                             % Komponentbalans över förstärkardelen
    i=i+1;
    [gamma1, gamma2] = wilson(x(i),W12,W21);
    Tstart = TB;
    TB=fsolve(@(T)find_Tb(T,x(i),gamma1,gamma2,Ant1,Ant2,P),Tstart);
    Psat1 = antoine(TB, Ant1);                                                % Ångtryck
    y(i) = (gamma1*x(i).*Psat1)/P;
end


%% del 2
% Värmen
Q_condensor = (Hvap*V*x(end)+Hvap*V*(1-x(end)))*10^3*3600^-1;      % W
Q_reboiler = (Hvap*v*y0 + Hvap*v*(1-y0))*10^3*3600^-1;             % W


% Faktorer
M_L = (x1*M1+(1-x1)*M2)*10^-3;      % kg mol^-1
M_V = (y0*M1+(1-y0)*M2)*10^-3;      % kg mol^-1

massandel_ethanol = (x1*46.07)/(x1*46.07+(1-x1)*60.0952);
rho_liquid = (massandel_ethanol*rho_ethanol+(1-massandel_ethanol)*rho_propanol);

F_LV = ((l*M_L)/(v*M_V))*sqrt(rho_vapour/rho_liquid);
C_F = 0.38*0.3048;                                         % Läs av diagram efter att F_LV beräknats (m/s)
F_F = 1;
F_ST = (liqsurften/20)^0.2;
F_HA = 1;
C = F_ST*F_F*F_HA*C_F;  

U_F = C*sqrt((rho_liquid-rho_vapour)/rho_vapour); 
U = 0.7*U_F; 

V_flowrate = v*10^3*3600^-1*M_V*rho_vapour^-1;
A_aktiv = V_flowrate/U; 
A = A_aktiv/0.8;
d = 2*sqrt(A/pi);


% UTSKRIVNING AV RESULTAT
disp([' ' ])
disp(['Ideala steg:               ' num2str(i)])
disp(['Diameter på tornet:        ' num2str(d)])
disp(['Höjd på tornet:            ' num2stri()])
disp([' ' ])

% Vapour and liquid flowrates in kmol h^-1
disp(['L (förstärkardel):         ' num2str(L)])
disp(['V (förstärkardel):         ' num2str(V)])
disp(['L (avdrivardel):           ' num2str(l)])
disp(['V (avdrivardell):          ' num2str(v)])
disp([' ' ])

% Värmen
disp(['Condenser duty (MW)        ' num2str(Q_condensor*10^-6)])
disp(['Reboiler duty (MW)         ' num2str(Q_reboiler*10^-6)])
disp([' ' ])

% Factors from correlation
disp(['F_LV                       ' num2str(F_LV)])
disp(['C_F (ft s^-1)              ' num2str(C_F/0.3048)])
disp([' ' ])


function [gamma1, gamma2] = wilson(x1,W12,W21)
    x2 = 1 - x1;

    gamma1 = exp(-log(x1+W12.*x2) + x2.*(W12./(x1+W12.*x2) - W21./(W21.*x1+x2)));
    gamma2= exp(-log(x2+W21.*x1) - x1.*(W12./(x1+W12.*x2) - W21./(W21.*x1+x2)));

end

function P_sat = antoine(T,Ant)

    A = Ant(1);
    B = Ant(2); 
    C = Ant(3);
    P_sat = exp(A-(B./(T+C))); 

end

function res = find_Tb(T,x1,gamma1,gamma2,Ant1,Ant2,P)

A1 = Ant1(1); B1 = Ant1(2); C1 = Ant1(3);
A2 = Ant2(1); B2 = Ant2(2); C2 = Ant2(3);

P1 = exp(A1-(B1./(T+C1)));
P2 = exp(A2-(B2./(T+C2)));

y1 = (x1.*gamma1.*P1)./P;
y2 = ((1-x1).*gamma2.*P2)./P;

res=y1+y2-1;

end