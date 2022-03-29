% Grundläggande kemiteknink, studioövning 5
% Värmeväxlare

clear, clc;

% Givna data
M_A=58.12*1000;  %[kg/mol]
M_B=56.11*1000;  %[kg/mol]
Fc_A=128*10^3;
                   %[mol/s]         Produktflöde (c=kalla flödet)
cpc=;       %[J/(molK)]      Produktflödets värmekapacitet
TCin=20;         %[C]            Produktflödets temperatur in
TCut=100;        %[C]            Produktflödets temperatur ut
mh=3.0;          %[mol/s]         Kondensatflöde (h=varma flödet)
cph=;      %[J/(molK)]      Kondensatflödets värmekapacitet
THin=90;         %[C]            Kondensatflödets temperatur in
U=1500;          %[W/(m2K)]      Värmegenomgångstal
Ka=600;          %[SEK/(m2 år)]  Kostnad för värmeväxlaren
beta=0.10e-3;    %[SEK/Wh]       Kostnad för ångan
tdrift=8760;     %[h/år]         Driftstid på ett år


% Beräkna kapacitetskoefficienterna, Cmin och Cmax.
C = [mc*cpc mh*cph];
Cmin = min(C);
Cmax = max(C);
%% 
%Beräkna optimala arean. (While-loop eller fsolve)
% Minns deriveringsregeln för en kvot: D (f/g) = (f'g-fg')/g^2

epsilon = @(NTU)(1 - exp(-NTU.*(1 - Cmin/Cmax)))/(1 - (Cmin/Cmax).*exp(-NTU.*(1 - Cmin/Cmax)));

NTU_guess = 7;
NTU = fsolve(@(NTU) ekonomi(NTU, Cmin, Cmax, Ka, U, THin, TCin, tdrift, beta), NTU_guess);

A = (NTU*Cmin)/U - Abef;

% Beräkna överfört värme i den nya värmeväxlaren.

qtot = epsilon(NTU).*Cmin.*(THin - TCin);
qbef = epsilon(U*Abef/Cmin).*Cmin.*(THin - TCin);
qnew = qtot - qbef;

%Beräkna andelen minskad ångeffekt.
T = qbef/(mc*cpc) + TCin;
Ts = qtot/(mc*cpc) + TCin;

qsteam_i = mc*cpc*(TCut - T);
qsteam_f = mc*cpc*(TCut - Ts);

steamsaved = qsteam_i - qsteam_f;

steamfrac = steamsaved/qsteam_i;



% Presentera resultaten.
disp(['Area nya VVX                     ' num2str(A, 3),' m^2'])
disp(['Sparad effekt av ånga            ' num2str(steamsaved*1e-3, 4),' kW'])
disp(['Q sparat relativt ångbehovet     ' num2str(steamfrac*100, 3),' %'])









