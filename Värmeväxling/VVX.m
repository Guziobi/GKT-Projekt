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
Tin=298;        %[K]            Produktfl�dets temperatur in
Tut_guess=750;  %[K]            Gissning av produktfl�dets temperatur ut f�r Cp-approx
Tmedel=(Tin+Tut_guess)./2;

%Produktfl�den
Fc_A=128*10^3;
Fc_B=5*10^3;
Fc_C=0;     
Fc_D=1091*10^3;  %[mol/h]

F_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molfl�de

%Molbr�k
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


%V�rmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);
cp_Dl=Cp_calc(Tmedel, CPcoeffDl);

%Total v�rmekapacitet f�r fl�det
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%Varma sidan
mvvx=2000;
Fvvx=1000*mvvx/M_D;       %[mol/h]
Tvvxin=800;             %[K]            Kondensatfl�dets temperatur in
cpvvx=Cp_calc(Tvvxin, CPcoeffD);   %[J/(kgK)]      Kondensatfl�dets v�rmekapacitet
         
U=1500;          %[W/(m2K)]      V�rmegenomg�ngstal

%Kostnader och �vrigt
Ka=600;          %[SEK/(m2 �r)]  Kostnad f�r v�rmev�xlaren
beta=0.10*10^-3; %[SEK/Wh]       Kostnad f�r �ngan
tdrift=8760;     %[h/�r]         Driftstid p� ett �r


%Ber�kning
C=[F_tot*cp_tot Fvvx*cpvvx];
Cmin=min(C);
Cmax=max(C);
Cmm=Cmin./Cmax;

%L�ser de/dNTU - bosse = 0, se ekonomi.m
NTU_guess = 7;
NTU = fzero(@(NTU) ekonomi(NTU, Cmin, Cmax, Cmm, Ka, U, Tvvxin, Tin, tdrift, beta), NTU_guess)

epsilon = (1 - exp(-NTU.*(1 - Cmin/Cmax)))/(1 - (Cmin/Cmax).*exp(-NTU.*(1 - Cmin/Cmax)));

A=NTU.*Cmin./U

Tcut=Tin+epsilon*Cmin*(Tvvxin-Tin)/(cp_tot*F_tot)


%% Kylning 5-6
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan f�r isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan f�r isobuten
M_D=18.01528;    %[kg/mol]       Molmassan f�r vatten

%Produktfl�de
%Temperaturer
Tin=600;        %[K]            Produktfl�dets temperatur in
Tut_guess=348;  %[K]            Gissning av produktfl�dets temperatur ut f�r Cp-approx
Tmedel=(Tin+Tut_guess)./2;

%Produktfl�den
Fc_A=20*10^3;
Fc_B=107*10^3;
Fc_C=107;     
Fc_D=1091*10^3;  %[mol/h]

F_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molfl�de

%Molbr�k
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


%V�rmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);
cp_Dl=Cp_calc(Tmedel, CPcoeffDl);

%Total v�rmekapacitet f�r fl�det
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%VVX-sidan
mvvx=2000;              %[kg/h]
Fvvx=1000*mvvx/M_D;     %[mol/h]
Tvvxin=287;             %[K]            Kondensatfl�dets temperatur in
cpvvx=Cp_calc(Tvvxin, CPcoeffDl);   %[J/(kgK)]      Kondensatfl�dets v�rmekapacitet
         
U=200;          %[W/(m2K)]      V�rmegenomg�ngstal, gas/v�tska-vvx

%Kostnader och �vrigt
a=32000; b=70; n=1.2; A_guess=1000;  %Kostnadsparametrar

Ka=a+b.*A_guess.^n.*9.99;   %[SEK/(m2 �r)]    Kostnad f�r v�rmev�xlaren (1 USD = 9.99 SEK)
beta=0.05;          %[SEK/Wh]       Kostnad f�r �ngan
tdrift=8000;        %[h/�r]         Driftstid p� ett �r

%Ber�kning
C=[F_tot*cp_tot Fvvx*cpvvx];
Cmin=min(C);
Cmax=max(C);
Cmm=Cmin./Cmax;

%L�ser de/dNTU - bosse = 0, se ekonomi.m
NTU_guess = 10;
NTU = fsolve(@(NTU) ekonomi(NTU, Cmin, Cmax, Cmm, Ka, U, Tvvxin, Tin, tdrift, beta), NTU_guess)

epsilon = 0.8;

A=NTU.*Cmin./U

Tut=Tin+epsilon*Cmin*(Tvvxin-Tin)/(cp_tot*F_tot)



%% Ugn 1-2
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan f�r isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan f�r isobuten
M_D=18.01528;    %[kg/mol]       Molmassan f�r vatten

%Produktfl�de
%Temperaturer
Tin=298;                %[K]            Produktfl�dets temperatur in
Tut=750;                %[K]            �nskad uttemperatur
Tmedel=(Tin+Tut)./2;    %[K]            Cp tas vid medeltemperatur        
dT=Tut-Tin;             %[K]            Temperaturdifferens

%Produktfl�den
Fc_A=128*10^3;
Fc_B=5*10^3;
Fc_C=0;     
Fc_D=1091*10^3;  %[mol/h]

F_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molfl�de

%Molbr�k
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

%V�rmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);
cp_Dl=Cp_calc(Tmedel, CPcoeffDl);

%Total v�rmekapacitet f�r produktfl�det
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%Produktfl�de
qh=cp_tot.*F_tot.*dT;    %[J/h]      V�rme per timme
q=qh/3600;               %[W]        V�rme per sekund

E=0.8;                   %Verkningsgrad ugn
qugn=(q/E)/10e6          %[MW] Ugnens effekt

%Kostnad
a=80000; b=109000; n=0.8;    %Kostnadsparametrar f�r cylindrisk ugn

K=a+b*qugn.^n*9.99          %[SEK/�r]


%% Ugn 3-4
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan f�r isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan f�r isobuten
M_D=18.01528;    %[kg/mol]       Molmassan f�r vatten

%Produktfl�de
%Temperaturer
Tin=500;                %[K]            Produktfl�dets temperatur in
Tut=750;                %[K]            �nskad uttemperatur
Tmedel=(Tin+Tut)./2;    %[K]            Cp tas vid medeltemperatur        
dT=Tut-Tin;             %[K]            Temperaturdifferens

%Produktfl�den
Fc_A=51*10^3;
Fc_B=77*10^3;
Fc_C=77*10^3;     
Fc_D=1091*10^3;  %[mol/h]

F_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/h]  Totalt molfl�de

%Molbr�k
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

%V�rmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);
cp_Dl=Cp_calc(Tmedel, CPcoeffDl);

%Total v�rmekapacitet f�r produktfl�det
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%Produktfl�de
qh=cp_tot.*F_tot.*dT;    %[J/h]      V�rme per timme
q=qh/3600;               %[W]        V�rme per sekund

E=0.8;                   %Verkningsgrad ugn
qugn=(q/E)/10e6          %[MW] Ugnens effekt

%Kostnad
a=80000; b=109000; n=0.8;    %Kostnadsparametrar f�r cylindrisk ugn

K=a+b*qugn.^n*9.99           %[SEK/�r]

%% Funktion f�r ber�kning av Cp(T)
function Cp = Cp_calc(T, CPcoeff)
    A = CPcoeff(1);
    B = CPcoeff(2);
    C = CPcoeff(3);
    D = CPcoeff(4);

    Cp = A + B*T + C*T.^2 + D*T.^3;

end