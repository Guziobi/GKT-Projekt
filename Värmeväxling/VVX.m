% Värmeväxling
% A=isobutan, B=isobuten, C=vätgas, D=vatten
%% VVX

clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan för isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan för isobuten
M_D=18.01528;    %[kg/mol]       Molmassan för vatten

%Kalla sidan (c)
%Temperaturer
Thin=298;        %[K]            Produktflödets temperatur in
Tut_guess=750;  %[K]            Gissning av produktflödets temperatur ut för Cp-approx
Tmedel=(Thin+Tut_guess)./2;

%Produktflöden
Fc_A=128*10^3;
Fc_B=5*10^3;
Fc_C=0;     
Fc_D=1091*10^3;  %[mol/h]

F_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molflöde

%Molbråk
yc_A=Fc_A./F_tot;
yc_B=Fc_B./F_tot;
yc_C=Fc_C./F_tot;
yc_D=Fc_D./F_tot;

% Cp-koefficienter
CPcoeffA = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];     % ISOBUTAN
CPcoeffB = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];     % ISOBUTEN
CPcoeffC = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];    % H2
CPcoeffD = [72.43 1.039*10^-2 -1.497*10^-6 0 ];         % H2O
CPcoeffDl = [32.24 0.001924 1.055e-5 -3.596e-9];        % H2O flytande


%Värmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);
cp_Dl=Cp_calc(Tmedel, CPcoeffDl);

%Total värmekapacitet för flödet
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%Varma sidan
mvvx=2000;
Fvvx=1000*mvvx/M_D;       %[mol/h]
Tcin=800;             %[K]            Kondensatflödets temperatur in
cpvvx=Cp_calc(Tcin, CPcoeffD);   %[J/(kgK)]      Kondensatflödets värmekapacitet
         
U=1500;          %[W/(m2K)]      Värmegenomgångstal

%Kostnader och övrigt
Ka=600;          %[SEK/(m2 år)]  Kostnad för värmeväxlaren
beta=0.10*10^-3; %[SEK/Wh]       Kostnad för ångan
tdrift=8760;     %[h/år]         Driftstid på ett år


%Beräkning
C=[F_tot*cp_tot Fvvx*cpvvx];
Cmin=min(C);
Cmax=max(C);
Cmm=Cmin./Cmax;

%Löser de/dNTU - bosse = 0, se ekonomi.m
NTU_guess = 7;
NTU = fzero(@(NTU) ekonomi(NTU, Cmin, Cmax, Cmm, Ka, U, Tcin, Thin, tdrift, beta), NTU_guess)

epsilon = (1 - exp(-NTU.*(1 - Cmin/Cmax)))/(1 - (Cmin/Cmax).*exp(-NTU.*(1 - Cmin/Cmax)));

A=NTU.*Cmin./U

Tcut=Thin+epsilon*Cmin*(Tcin-Thin)/(cp_tot*F_tot)


%% Kylning 5-6
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan för isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan för isobuten
M_D=18.01528;    %[kg/mol]       Molmassan för vatten

%Produktflöde
%Temperaturer
Thin=600;        %[K]            Produktflödets temperatur in
Tut_guess=348;  %[K]            Önskad uttemp för Cp-approx, manuell iterering till denna temp nås
Tmedel=(Thin+Tut_guess)./2;    %[K] Cp tas vid medeltemp

%Produktflöden [mol/s]
Fc_A=20*10^3./3600;
Fc_B=107*10^3./3600;
Fc_C=107./3600;     
Fc_D=1091*10^3./3600;

F_tot=(Fc_A+Fc_B+Fc_C+Fc_D);  %[mol/s]  Totalt molflöde

%Molbråk
yc_A=Fc_A./F_tot;
yc_B=Fc_B./F_tot;
yc_C=Fc_C./F_tot;
yc_D=Fc_D./F_tot;

% Cp-koefficienter
CPcoeffA = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];     % ISOBUTAN
CPcoeffB = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];     % ISOBUTEN
CPcoeffC = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];    % H2
CPcoeffD = [72.43 1.039*10^-2 -1.497*10^-6 0 ];         % H2O
CPcoeffDl = [32.24 0.001924 1.055e-5 -3.596e-9];        % H2O flytande

%Värmekapaciteter [J/(mol K)]
cp_A=Cp_calc(Tmedel, CPcoeffA);        
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);
cp_Dl=Cp_calc(Tmedel, CPcoeffDl);

%Total värmekapacitet för flödet
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%VVX-sidan (kalla)
mvvx=54000/3600;              %[kg/s] Massflöde
Fvvx=(1000*mvvx/M_D);     %[mol/s] Molflöde
Tcin=287;             %[K] Kondensatflödets temperatur in
cpvvx=Cp_calc(Tcin, CPcoeffDl);   %[J/(mol K)]      Kondensatflödets värmekapacitet
        
U=200;                  %[W/(m2K)]      Värmegenomgångstal, gas/vätska-vvx

%Kostnader och övrigt
a=32000; b=70; n=1.2;  %Kostnadsparametrar

A_guess=492;          %Gissning av area för iterering
                       %Arean bestämmer kostnaden som i sin tur används för
                       %att optimera arean

Ka=a+b.*A_guess.^n.*9.99;   %[SEK/(m2 år)]    Kostnad för värmeväxlaren (1 USD = 9.99 SEK)
beta=0.05;                  %[SEK/Wh]       Kostnad för ångan
tdrift=8000;                %[h/år]         Driftstid på ett år

%Beräkning
C=[F_tot*cp_tot Fvvx*cpvvx];
Cmin=min(C);
Cmax=max(C);
Cmm=Cmin./Cmax;         %Kvoten Cmin/Cmax

%Optimering av area likt studioövning 5
%Löser de/dNTU - bosse = 0, se funktionsfilen ekonomi.m
NTU_guess = 4;
NTU = fsolve(@(NTU) ekonomi(NTU, Cmm, Ka, U, Tcin, Thin, tdrift, beta), NTU_guess)


epsilon = 0.8;

A=NTU.*Cmin./U

%Uttemp beräknad med värmebalans
Tcut=Thin+epsilon*Cmin*(Tcin-Thin)/(cp_tot*F_tot)


%% Ugn 1-2
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan för isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan för isobuten
M_D=18.01528;    %[kg/mol]       Molmassan för vatten

%Produktflöde
%Temperaturer
Thin=298;                %[K]            Produktflödets temperatur in
Tut=750;                %[K]            Önskad uttemperatur
Tmedel=(Thin+Tut)./2;    %[K]            Cp tas vid medeltemperatur        
dT=Tut-Thin;             %[K]            Temperaturdifferens

%Produktflöden
Fc_A=128*10^3;
Fc_B=5*10^3;
Fc_C=0;     
Fc_D=1091*10^3;  %[mol/h]

F_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molflöde

%Molbråk
yc_A=Fc_A./F_tot;
yc_B=Fc_B./F_tot;
yc_C=Fc_C./F_tot;
yc_D=Fc_D./F_tot;

% Cp-koefficienter
CPcoeffA = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];     % ISOBUTAN
CPcoeffB = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];     % ISOBUTEN
CPcoeffC = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];    % H2
CPcoeffD = [72.43 1.039*10^-2 -1.497*10^-6 0 ];         % H2O
CPcoeffDl = [32.24 0.001924 1.055e-5 -3.596e-9];        % H2O flytande

%Värmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);
cp_Dl=Cp_calc(Tmedel, CPcoeffDl);

%Total värmekapacitet för produktflödet
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%Produktflöde
qh=cp_tot.*F_tot.*dT;    %[J/h]      Värme per timme
q=qh/3600;               %[W]        Värme per sekund

E=0.8;                   %Verkningsgrad ugn
qugn=(q/E)/10e6          %[MW] Ugnens effekt

%Kostnad
a=80000; b=109000; n=0.8;    %Kostnadsparametrar för cylindrisk ugn

K=a+b*qugn.^n*9.99          %[SEK/år]


%% Ugn 3-4
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan för isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan för isobuten
M_D=18.01528;    %[kg/mol]       Molmassan för vatten

%Produktflöde
%Temperaturer
Thin=500;                %[K]            Produktflödets temperatur in
Tut=750;                %[K]            Önskad uttemperatur
Tmedel=(Thin+Tut)./2;    %[K]            Cp tas vid medeltemperatur        
dT=Tut-Thin;             %[K]            Temperaturdifferens

%Produktflöden
Fc_A=51*10^3;
Fc_B=77*10^3;
Fc_C=77*10^3;     
Fc_D=1091*10^3;  %[mol/h]

F_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molflöde

%Molbråk
yc_A=Fc_A./F_tot;
yc_B=Fc_B./F_tot;
yc_C=Fc_C./F_tot;
yc_D=Fc_D./F_tot;

% Cp-koefficienter
CPcoeffA = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];     % ISOBUTAN
CPcoeffB = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];     % ISOBUTEN
CPcoeffC = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];    % H2
CPcoeffD = [72.43 1.039*10^-2 -1.497*10^-6 0 ];         % H2O
CPcoeffDl = [32.24 0.001924 1.055e-5 -3.596e-9];        % H2O flytande

%Värmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);
cp_Dl=Cp_calc(Tmedel, CPcoeffDl);

%Total värmekapacitet för produktflödet
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%Produktflöde
qh=cp_tot.*F_tot.*dT;    %[J/h]      Värme per timme
q=qh/3600;               %[W]        Värme per sekund

E=0.8;                   %Verkningsgrad ugn
qugn=(q/E)/10e6          %[MW] Ugnens effekt

%Kostnad
a=80000; b=109000; n=0.8;    %Kostnadsparametrar för cylindrisk ugn

K=a+b*qugn.^n*9.99           %[SEK/år]

%% Funktion för beräkning av Cp(T)
function Cp = Cp_calc(T, CPcoeff)
    A = CPcoeff(1);
    B = CPcoeff(2);
    C = CPcoeff(3);
    D = CPcoeff(4);

    Cp = A + B*T + C*T.^2 + D*T.^3;

end