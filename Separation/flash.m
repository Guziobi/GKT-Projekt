%Flashtankar
clc, clear

%todo: Kika enheter!

%Data
k = 0.107; %m s-1 @8bar
% Ångflöde in
V = 1300; %kmol/h
%Vätskeflöde från tanken
L = 1100; %kmol/h
%uppehållstid
tau = 1/6; % h (10min)

rho_V = 554; %placeholder
rho_L = 559; %placeholder

ut = k*sqrt((rho_L - rho_V)/rho_V); %m s-1

% Diameter
D = sqrt(4*V/(pi*0.15*ut)); %m

% Vätskehöjd
HL = (L*tau)/((pi*D^2)/4);

%Höjd flashtank
H = HL + 1.5*D;