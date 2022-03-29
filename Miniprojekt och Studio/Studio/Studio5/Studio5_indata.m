% Grundl�ggande kemiteknink, studio�vning 5
% V�rmev�xlare

clear, clc;

% Givna data
M_A=58.12*1000;  %[kg/mol]
M_B=56.11*1000;  %[kg/mol]
Fc_A=128*10^3;
                   %[mol/s]         Produktfl�de (c=kalla fl�det)
cpc=;       %[J/(molK)]      Produktfl�dets v�rmekapacitet
TCin=20;         %[C]            Produktfl�dets temperatur in
TCut=100;        %[C]            Produktfl�dets temperatur ut
mh=3.0;          %[mol/s]         Kondensatfl�de (h=varma fl�det)
cph=;      %[J/(molK)]      Kondensatfl�dets v�rmekapacitet
THin=90;         %[C]            Kondensatfl�dets temperatur in
U=1500;          %[W/(m2K)]      V�rmegenomg�ngstal
Ka=600;          %[SEK/(m2 �r)]  Kostnad f�r v�rmev�xlaren
beta=0.10e-3;    %[SEK/Wh]       Kostnad f�r �ngan
tdrift=8760;     %[h/�r]         Driftstid p� ett �r


% Ber�kna kapacitetskoefficienterna, Cmin och Cmax.
C = [mc*cpc mh*cph];
Cmin = min(C);
Cmax = max(C);
%% 
%Ber�kna optimala arean. (While-loop eller fsolve)
% Minns deriveringsregeln f�r en kvot: D (f/g) = (f'g-fg')/g^2

epsilon = @(NTU)(1 - exp(-NTU.*(1 - Cmin/Cmax)))/(1 - (Cmin/Cmax).*exp(-NTU.*(1 - Cmin/Cmax)));

NTU_guess = 7;
NTU = fsolve(@(NTU) ekonomi(NTU, Cmin, Cmax, Ka, U, THin, TCin, tdrift, beta), NTU_guess);

A = (NTU*Cmin)/U - Abef;

% Ber�kna �verf�rt v�rme i den nya v�rmev�xlaren.

qtot = epsilon(NTU).*Cmin.*(THin - TCin);
qbef = epsilon(U*Abef/Cmin).*Cmin.*(THin - TCin);
qnew = qtot - qbef;

%Ber�kna andelen minskad �ngeffekt.
T = qbef/(mc*cpc) + TCin;
Ts = qtot/(mc*cpc) + TCin;

qsteam_i = mc*cpc*(TCut - T);
qsteam_f = mc*cpc*(TCut - Ts);

steamsaved = qsteam_i - qsteam_f;

steamfrac = steamsaved/qsteam_i;



% Presentera resultaten.
disp(['Area nya VVX                     ' num2str(A, 3),' m^2'])
disp(['Sparad effekt av �nga            ' num2str(steamsaved*1e-3, 4),' kW'])
disp(['Q sparat relativt �ngbehovet     ' num2str(steamfrac*100, 3),' %'])









