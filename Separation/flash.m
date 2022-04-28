%Flashtankar
clc, clear

%todo: Kika enheter!

%Data
k = 0.107; %m s-1 @8bar
% Kokvarm ström in
T = 1558.5; %kmol/h
%Vätskeflöde från tanken
L = 1310; %kmol/h
% Ångflöde ut
V = T-L; %kmol/h

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


%% tank 2
T2 = V;
V2 = 32.6381;
L2 = T2-V2;

disp(['Diameter tank 1:               ' num2str(D)])
disp(['Diameter på tornet:        ' num2str(d)])
disp(['Höjd på tornet:            ' num2str(h)])
disp([' ' ])
