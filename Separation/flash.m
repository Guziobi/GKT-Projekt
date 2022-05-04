%Flashtankar

%todo: Kika enheter!

% Data
k = 0.107; %m s-1 @8bar
p = 101325; % Pa
R = 8.3145;
T = 348;

%Antoinekonstanter A  B  C
Ant1 =  [15.7564 2132.42 -33.15];  % buten
Ant2 =  [15.6782 2154.90 -34.42];  % butan 
%Wilsonfaktorer
W12 = 0.48584; 
W21 = 1.64637; 

% Molmassor
MA = 56.1063;      % g mol-1  Butan
MB = 58.1222;      % g mol-1  Buten
MC = 2.0160;       % g mol-1  Vätgas
MD = 18.0160;      % g mol-1  Vatten

rhoA = 559.0;       % kgm-3 Buten
rhoB = 556.62;      % kgm-3 Butan
rhoD = 974.93;      % kg m-3 Vätska

%V =nRT/p
% Vätskeflöde från tanken (Vatten)
L = (U(end,4)*MD*1e-3)/rhoD; %m3/s
% Ångflöde ut
V = ((sum(U(end,:))-U(end,4)-U(end,5))*R*T)/p; %m3/s

%uppehållstid
tau = 10*60; % s (10min)
rho_L = rhoD;

%sum(xi*(MiP/RT))
rho_V = p*1e-3*(U(end,1)*MA + U(end,2)*MB + U(end,3)*MC)/(T*R*V);
ut = k*sqrt((rho_L - rho_V)/rho_V); %m s-1

% Diameter
D1 = sqrt(4*V/(pi*0.15*ut)); %m

% Vätskehöjd
HL1 = (L*tau)/((pi*D1.^2)/4);

%Höjd flashtank
H1 = HL1 + 1.5*D1;


%% tank 2
p = 101325*6; % Pa

x1 = U(end,2)/(U(end,1) + U(end,2));

[gamma1, gamma2] = wilson(x1,W12,W21);
Tguess = 273;
T = fsolve(@(T) find_Tb(T,x1,gamma1,gamma2,Ant1,Ant2,p*0.0075006168), Tguess);

V2 = (U(end,3)*R*T)/p; %m3/s
L2 = (U(end,1)*MA*1e-3)/rhoA + (U(end,2)*MB*1e-3)/rhoB; %m3/s

rho_L = (U(end,1)*rhoA + U(end,2)*rhoB)/(U(end,1) + U(end,2)); %kg m-3
rho_V = p*1e-3*(U(end,3)*MC)/(T*R*V2); %kg m-3

ut = k*sqrt((rho_L - rho_V)/rho_V); %m s-1

% Diameter
D2 = sqrt(4*V2/(pi*0.15*ut)); %m

% Vätskehöjd
HL2 = (L2*tau)/((pi*D2.^2)/4);

%Höjd flashtank
H2 = HL2 + 1.5*D2;


disp(['Diameter tank 1:               ' num2str(D1)])
disp(['Höjd på tank 1:            ' num2str(H1)])
disp([' ' ])

disp(['Diameter tank 2:               ' num2str(D2)])
disp(['Höjd på tank 2:            ' num2str(H2)])
disp([' ' ])

function res = find_Tb(T,x1,gamma1,gamma2,Ant1,Ant2,P)

A1 = Ant1(1); B1 = Ant1(2); C1 = Ant1(3);
A2 = Ant2(1); B2 = Ant2(2); C2 = Ant2(3);

P1 = exp(A1-(B1./(T+C1)));
P2 = exp(A2-(B2./(T+C2)));

y1 = (x1.*gamma1.*P1)./P;
y2 = ((1-x1).*gamma2.*P2)./P;

res=y1+y2-1;

end

function [gamma1, gamma2] = wilson(x1,W12,W21)
    x2 = 1 - x1;

    gamma1 = exp(-log(x1+W12.*x2) + x2.*(W12./(x1+W12.*x2) - W21./(W21.*x1+x2)));
    gamma2= exp(-log(x2+W21.*x1) - x1.*(W12./(x1+W12.*x2) - W21./(W21.*x1+x2)));

end





