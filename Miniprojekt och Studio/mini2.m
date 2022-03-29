clc, clear
% Indata
T0 = 502; % K
Tref = 273.15; % K
q0 = 21e-6; % m^3 s^-1
Ptot = 1.2e5; % Pa
Ea = 115e3; % J mol^-1
A = 3e4; % m^3 mol^-1 s^-1
DHr = -170e3; %J mol^-1
Xr = 0.4; % 40% conversion
% Cp-värden J mol^-1 K^-1
Cp_c2h4 = 80;
Cp_c4h6 = 150;
Cp_c6h10 = 230;
Cp_steam = 36;

%Övriga konstanter
gaskonstant = 8.31446261815324; % m^3 Pa K^-1 mol^-1


%% Different flows mol/s
F0_tot = (Ptot*q0)/(gaskonstant*T0);
F0_c2h4 = 0.25*F0_tot;
F0_c4h6 = F0_c2h4;
F0_steam = 0.5*F0_tot;

R = @(X) X.*F0_c2h4;
F1_steam = F0_steam;
F1_c2h4 = @(X) F0_c2h4 - X.*F0_c2h4;
F1_c4h6 = @(X) F0_c4h6 - X.*F0_c4h6;
F1_c6h10 = @(X) X.*F0_c2h4;
F1_tot = @(X) F1_steam + F1_c2h4(X) + F1_c4h6(X) + F1_c6h10(X);

%% Total heat effect of reaction in W (J/s)
Qreak = @(X) F0_c2h4.*X.*(-DHr);

%% Mera värme
sum_F0_Cp = F0_c2h4.*Cp_c2h4 + F0_c4h6.*Cp_c4h6 + F0_steam.*Cp_steam;
sum_F1_Cp = @(X) F1_c2h4(X).*Cp_c2h4 + F1_c4h6(X).*Cp_c4h6 + F1_steam(X).*Cp_steam + F1_c6h10(X).*Cp_c6h10;

Tut = @(X) Qreak(X)./sum_F0_Cp + T0;

%% Reaktion rate ideal tank

% k1 från Arrhenius
k = @(X) A.*exp(-Ea./(gaskonstant.*Tut(X))); % m^3 mol^-1 s^-1
k(0.4)

% C1_c2h4 = @(X) (Ptot./(gaskonstant.*Tut(X))).*((1-X)./(4-X)); % mol dm^-3
C1_c2h4 = @(X) (Ptot./(gaskonstant.*Tut(X))).*((F0_c4h6.*(1-X))./(F1_tot(X))); % mol dm^-3

r1 = @(X) k(X).*C1_c2h4(X).^2;

%% Volym CSTR
V_cstr = @(X) (X*F0_c2h4)./r1(X);

%% Värden intressanta värden till mobius

F0_c2h4*1e3 

Qreak = Qreak(Xr)

T_ut = Tut(Xr)

r = r1(Xr)*1e3

V_cstr = V_cstr(Xr)

%% Volym PFR
V_PFR = F0_c2h4*integral(@(X) 1./r1(X), 0, Xr)













