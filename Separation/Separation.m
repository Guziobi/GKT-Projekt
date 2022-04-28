% Separation
clc, clear
options = optimset('Display','off');    % Så att skit inte skrivs ut efter fsolve

% Data - separation 
q = 1;             % kokvarmt tillflöde
F = 131;           % kmol h-1
P = 2280;          % mmHg (3 atm)
xf = 0.8969;         % molbråk buten
xd = 0.95;         % destillatbråk 
xb = 0.05;         % bottenbråk
R = 5;            % återflödesförhållande
%Molmassor
M1 = 56.1063;      % g mol-1
M2 = 58.1222;      % g mol-1
% Molmassor
M1 = 56.1063;       % g mol-1
M2 = 58.1222;       % g mol-1
% Kokpunkter
Tb_1 = -6.3+273.15; % K   (buten)
Tb_2 = -1+273.15;   % K   (butan)


%Antoinekonstanter A  B  C
Ant1 =  [15.7564 2132.42 -33.15];  % buten
Ant2 =  [15.6782 2154.90 -34.42];  % butan 
%Wilsonfaktorer
W12 = 0.48584; 
W21 = 1.64637; 
%Densiteter
L_rho1 = 559.0;          % kgm-3 Buten
L_rho2 = 556.62;         % kgm-3 Butan
%Ytspänning
surfaceten = 24;          % dyn cm-1
    
% Beräkning av flöden 
B = F*((xf-xd)/(xb-xd));
D = F-B;

L = R*D; 
V = L+D;

l = L + q*F;    % L-streck
v = l-B;        % V-streck

% Återkokare 
T = zeros(1,30);
[gamma1, gamma2] = wilson(xb,W12,W21);
Tstart = 273.15;
TB=fsolve(@(T)find_Tb(T,xb,gamma1,gamma2,Ant1,Ant2,P),Tstart,options);
TB_reboiler = TB;
Psat1 = antoine(TB, Ant1);
y0 = (gamma1*xb.*Psat1)/P;
x1 = (v/l)*y0 + (L_rho1/l)*xb;

% Avdrivardel 
x = zeros(1,60);
y = zeros(1,60);
x(1) = x1;
y(1)=y0;
i = 0; 

% Jämviktkurva och jämviktsplot
xeq = 0:0.001:1;    
[gamma1, gamma2] = wilson(xeq,W12,W21);
Tstart = linspace(-6.3+273.15,-1+273.15,1001);
TBeq=fsolve(@(T)find_Tb(T,xeq,gamma1,gamma2,Ant1,Ant2,P),Tstart,options);
Psat1 = antoine(TBeq, Ant1);                                            % Ångtryck
yeq = (gamma1.*xeq.*Psat1)./P;

plot(xeq,xeq)   % Referenslinje
hold on
plot(xeq,yeq)
xlabel('x_1'), ylabel('y_1')
axis([0 1 0 1])
legend('Jämviktskurva', 'Referenslinje','Location','northwest')

% Sorels metod
while x<xf
    i = i+1; 
    [gamma1, gamma2] = wilson(x(i),W12,W21);
    Tstart = TB;
    TB=fsolve(@(T)find_Tb(T,x(i),gamma1,gamma2,Ant1,Ant2,P),Tstart,options);
    T(i) = TB;
    Psat1 = antoine(TB, Ant1);                                         % Ångtryck
    y(i) = (gamma1*x(i).*Psat1)/P;
    x(i+1) = (v/l)*y(i) + (B/l)*xb;
end

m = i+1;

% Förstärkare
while y(i)<xd
    x(i+1)= (V/L)*y(i) + (1/L)*(B*xb-F*xf);                        % Komponentbalans över förstärkardelen
    i=i+1;
    [gamma1, gamma2] = wilson(x(i),W12,W21);
    Tstart = TB;
    TB=fsolve(@(T)find_Tb(T,x(i),gamma1,gamma2,Ant1,Ant2,P),Tstart,options);
    T(i) = TB;
    Psat1 = antoine(TB, Ant1);                                      % Ångtryck
    y(i) = (gamma1*x(i).*Psat1)/P;
end

bottnar = i/0.7;

%% Dimensionering

% Sammansättning ut ur återkokare
v_x1 = y0;              % Buten
v_x2 = 1-y0;            % Butan
% Sammansättning in i återkokare
l_x1 = x(1);
l_x2 = 1 - x(1);

% Flödesfaktorer
surftention = 70; %dyn cm-1
Fst = (surftention/20)^0.2;
vaporveloc = 0.7;
Ff = 1; %nono-foaming
Fha = 1; %Hålen är bra

%Bottenavstånd
trayheight = 0.45; %m

%Densiteter för vätska och gas
rho_L = ((l_x1*M1*1e-3)/(l_x1*M1*1e-3 + l_x2*M2*1e-3))*L_rho1 + ((l_x2*M2*1e-3)/(l_x1*M1*1e-3 + l_x2*M2*1e-3))*L_rho2; % kg m-3
rho_V = v_x1*M1*1e-3*(P*133.322368/(R*TB_reboiler)) + v_x2*M2*1e-3*(P*133.322368/(R*TB_reboiler));

%molmassor
M_L = l_x1*M1 + l_x2*M2;
M_V = v_x1*M1 + v_x2*M2;

%Belastningsparameter
Flv = ((l*M_L)/(v*M_V)) * sqrt(rho_V/rho_L);

% från diagramet
Cf = 0.25*0.3048; %ft/s -> m/s

%Flödningsparametern
C = Fst*Ff*Fha*Cf;

% Ånghastigheten vid flödning
Uf = C*sqrt((rho_L-rho_V)/rho_V);

%Ånghastighet
U = vaporveloc*Uf;

%Aktiv area
Aaktiv = (V*(1/3.6)*M_V*1e-3*(1/rho_V))/U;

%Total area
Atot = Aaktiv/0.8;

%Bottendiameter
d = 2*sqrt(Atot/pi);

%Kolonnens höjd
h = trayheight * (bottnar + 1);

%% Värmen

Hvap1 = 20.6e3;  % j mol-1
Hvap2 = 19.99e3; % j mol-1

% Värmen
Q_condensor = (Hvap1*V*1e3*x(end) + Hvap2*V*1e3*(1-x(end))) * 3600^-1;    % W
Q_reboiler = (Hvap1*v*1e3*y0 + Hvap2*v*1e3*(1-y0)) * 3600^-1;             % W

%% Kostnader

% Parametrar givna i PM
Param = [130 440 1.8        % sieve tray
         210 400 1.9        % valve tray
         340 400 1.9];      % bubble cap tray
     
kurs = 9.99;                % Växelkursen sek/dollar
lang = 4;                   % Langfaktorn
index = 596/532.9;
     
kostnad_sieve = Cost(d,Param(1,:))*kurs*lang*index*bottnar
kostnad_valve = Cost(d,Param(2,:))*kurs*lang*index*bottnar
kostnad_bubble = Cost(d,Param(3,:))*kurs*lang*index*bottnar


%% UTSKRIVNING AV RESULTAT
disp([' ' ])
disp(['Ideala steg:               ' num2str(bottnar)])
disp(['Diameter på tornet:        ' num2str(d)])
disp(['Höjd på tornet:            ' num2str(h)])
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
disp(['F_LV                       ' num2str(Flv)])
disp(['C_F (ft s^-1)              ' num2str(Cf/0.3048)])
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

function kostnad = Cost(S,Param)
    a = Param(1);
    b = Param(2);
    n = Param(3);
    
    kostnad = a + b*S.^n;
end