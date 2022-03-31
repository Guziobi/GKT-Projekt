% V�rmev�xling
% A=isobutan, B=isobuten, C=v�tgas, D=vatten
%% VVX 1

clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan f�r isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan f�r isobuten

%Kalla sidan (c)
%Temperaturer
Tcin=298;        %[K]            Produktfl�dets temperatur in?
Tcut=750;        %[K]            Produktfl�dets temperatur ut?

%Produktfl�den
Fc_A=128*10^3;
Fc_B=5*10^3;
Fc_C=0;     
Fc_D=1091*10^3;  %[mol/h]        Produktfl�de (c=kalla fl�det)??

Fc_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molfl�de

%Molbr�k
yc_A=Fc_A./Fc_tot;
yc_B=Fc_B./Fc_tot;
yc_C=Fc_C./Fc_tot;
yc_D=Fc_D./Fc_tot;

%V�rmekapaciteter
cpc_A=cp_A(Tcin);
cpc_B=cp_B(Tcin);
cpc_C=cp_C(Tcin);
cpc_D=cp_D(Tcin);

cpc_tot=yc_A.*cpc_A+yc_B.*cpc_B+yc_C.*cpc_C+yc_D.*cpc_D; %[J/(molK)]  Total v�rmekapacitet

%Varma sidan
cph=4.18*10^3;   %[J/(kgK)]      Kondensatfl�dets v�rmekapacitet
Thin=90;         %[C]            Kondensatfl�dets temperatur in
U=1500;          %[W/(m2K)]      V�rmegenomg�ngstal

%Kostnader och �vrigt
Ka=600;          %[SEK/(m2 �r)]  Kostnad f�r v�rmev�xlaren
beta=0.10*10^-3; %[SEK/Wh]       Kostnad f�r �ngan
tdrift=8760;     %[h/�r]         Driftstid p� ett �r


epsilon = @(NTU)(1 - exp(-NTU.*(1 - Cmin/Cmax)))/(1 - (Cmin/Cmax).*exp(-NTU.*(1 - Cmin/Cmax)));



