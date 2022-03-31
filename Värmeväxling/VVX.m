% Värmeväxling
% A=isobutan, B=isobuten, C=vätgas, D=vatten
%% VVX 1

clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan för isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan för isobuten

%Kalla sidan (c)
%Temperaturer
Tcin=298;        %[K]            Produktflödets temperatur in?
Tcut=750;        %[K]            Produktflödets temperatur ut?
Tcmedel=(Tcin+Tcut)./2;

%Produktflöden
Fc_A=128*10^3;
Fc_B=5*10^3;
Fc_C=0;     
Fc_D=1091*10^3;  %[mol/h]        Produktflöde (c=kalla flödet)??

Fc_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molflöde

%Molbråk
yc_A=Fc_A./Fc_tot;
yc_B=Fc_B./Fc_tot;
yc_C=Fc_C./Fc_tot;
yc_D=Fc_D./Fc_tot;

%Värmekapaciteter
cpc_A=cp_A(Tcmedel);
cpc_B=cp_B(Tcmedel);
cpc_C=cp_C(Tcmedel);
cpc_D=cp_D(Tcmedel);

cpc_tot=yc_A.*cpc_A+yc_B.*cpc_B+yc_C.*cpc_C+yc_D.*cpc_D; %[J/(molK)]  Total värmekapacitet för flödet

%Varma sidan
Fh=10000;        %[mol/h]
Thin=800;        %[K]            Kondensatflödets temperatur in
cph=cp_Dl(Thin); %[J/(kgK)]      Kondensatflödets värmekapacitet
         
U=1500;          %[W/(m2K)]      Värmegenomgångstal

%Kostnader och övrigt
Ka=600;          %[SEK/(m2 år)]  Kostnad för värmeväxlaren
beta=0.10*10^-3; %[SEK/Wh]       Kostnad för ångan
tdrift=8760;     %[h/år]         Driftstid på ett år


%Beräkning
C=[Fc_tot*cpc_tot Fh*cph];
Cmin=min(C);
Cmax=max(C);

epsilon=Fc_tot.*cpc_tot.*(Tcut-Tcin)./(Cmin.*(Thin-Tcin))
%epsilon = @(NTU)(1 - exp(-NTU.*(1 - Cmin/Cmax)))/(1 - (Cmin/Cmax).*exp(-NTU.*(1 - Cmin/Cmax)));




