% V�rmev�xling
% A=isobutan, B=isobuten, C=v�tgas, D=vatten
%% VVX 1

clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan f�r isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan f�r isobuten
M_D=18.01528;    %[kg/mol]       Molmassan f�r vatten

%Kalla sidan (c)
%Temperaturer
Tcin=298;        %[K]            Produktfl�dets temperatur in
Tcut_guess=750;  %[K]            Gissning av produktfl�dets temperatur ut
Tcmedel=(Tcin+Tcut_guess)./2;

%Produktfl�den
Fc_A=128*10^3;
Fc_B=5*10^3;
Fc_C=0;     
Fc_D=1091*10^3;  %[mol/h] 

Fc_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molfl�de

%Molbr�k
yc_A=Fc_A./Fc_tot;
yc_B=Fc_B./Fc_tot;
yc_C=Fc_C./Fc_tot;
yc_D=Fc_D./Fc_tot;

%V�rmekapaciteter
cpc_A=cp_A(Tcmedel);        %[J/(molK)]
cpc_B=cp_B(Tcmedel);
cpc_C=cp_C(Tcmedel);
cpc_D=cp_D(Tcmedel);

cpc_tot=yc_A.*cpc_A+yc_B.*cpc_B+yc_C.*cpc_C+yc_D.*cpc_D; %[J/(molK)]  Total v�rmekapacitet f�r fl�det

%Varma sidan
mh=1800;
Fh=1000*mh/M_D;       %[mol/h]
Thin=800;        %[K]            Kondensatfl�dets temperatur in
cph=cp_Dl(Thin); %[J/(kgK)]      Kondensatfl�dets v�rmekapacitet
         
U=1500;          %[W/(m2K)]      V�rmegenomg�ngstal

%Kostnader och �vrigt
Ka=600;          %[SEK/(m2 �r)]  Kostnad f�r v�rmev�xlaren
beta=0.10*10^-3; %[SEK/Wh]       Kostnad f�r �ngan
tdrift=8760;     %[h/�r]         Driftstid p� ett �r


%Ber�kning
C=[Fc_tot*cpc_tot Fh*cph];
Cmin=min(C);
Cmax=max(C);
Cmm=Cmin./Cmax;

%L�ser de/dNTU - bosse = 0, se ekonomi.m
NTU_guess = 7;
NTU = fzero(@(NTU) ekonomi(NTU, Cmin, Cmax, Cmm, Ka, U, Thin, Tcin, tdrift, beta), NTU_guess)

epsilon = (1 - exp(-NTU.*(1 - Cmin/Cmax)))/(1 - (Cmin/Cmax).*exp(-NTU.*(1 - Cmin/Cmax)));

A=NTU.*Cmin./U

Tcut=Tcin+epsilon*Cmin*(Thin-Tcin)/(cpc_tot*Fc_tot)
