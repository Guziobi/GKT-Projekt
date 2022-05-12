%% Förvärming med produktflöde
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan för isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan för isobuten
M_D=18.01528;    %[kg/mol]       Molmassan för vatten

%Temperaturer
Thin=929;            %[K]            Produktflödets intemperatur (varm sida)
Tcin=180+273;        %[K]            Reaktantflödets intemperatur (kall sida) (blandat med 180C ånga)
Tmedel=(Thin+Tcin)./2;     %[K] Cp tas vid medeltemp

%Flöde varm sida [mol/s]
Fh_A=4.4146;
Fh_B=32.5298;
Fh_C=32.5298;     
Fh_D=1900*10^3./3600;

Fh_tot=(Fh_A+Fh_B+Fh_C+Fh_D);  %[mol/s]  Totalt molflöde varm sida

%Flöden kall sida
Fc_A=190*10^3./3600;
Fc_B=0;
Fc_C=0;     
Fc_D=1900*10^3./3600;

Fc_tot=(Fc_A+Fc_B+Fc_C+Fc_D);  %[mol/s]  Totalt molflöde varm sida

%Molbråk varm sida
yh_A=Fh_A./Fh_tot;
yh_B=Fh_B./Fh_tot;
yh_C=Fh_C./Fh_tot;
yh_D=Fh_D./Fh_tot;

%Molbråk kall sida
yc_A=Fc_A./Fc_tot;
yc_B=Fc_B./Fc_tot;
yc_C=Fc_C./Fc_tot;
yc_D=Fc_D./Fc_tot;

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

%Total värmekapacitet, varm sida
cph_tot=yh_A.*cp_A+yh_B.*cp_B+yh_C.*cp_C+yh_D.*cp_D; %[J/(molK)] 

%Total värmekapacitet, varm sida
cpc_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

U=50;                   %[W/(m2 K)]      Värmegenomgångstal, gas/gas-vvx
A=800;                  %[m2] Area
%A=linspace(0,1000,1000)    %För att plotta epsilon mot A

%Kostnader
a=32000; b=70; n=1.2;   %Kostnadsparametrar
Ka=a+b.*A.^n.*9.99;     %[SEK]    Kostnad för värmeväxlaren (1 USD = 9.99 SEK)

%Beräkning
C=[Fh_tot*cph_tot Fc_tot*cpc_tot];
Cmin=min(C);
Cmax=max(C);
Cmm=Cmin./Cmax;         %Kvoten Cmin/Cmax

NTU=U*A/Cmin;           %Number of transfered units


epsilon=(1-exp(-NTU.*(1-Cmm)))./(1-Cmm.*exp(-NTU.*(1-Cmm)));    %Termisk verkningsgrad


%Sökt
q=epsilon*Cmin*(Thin-Tcin);          %[W] Effekt
Tcut=Tcin+q./(Fc_tot.*cpc_tot);      %[K] Uttemperatur kall sida
Thut=Thin-q./(Fh_tot.*cph_tot);      %[K] Uttemperatur varm sida

%plot(A,epsilon) %För optimering av area

disp(['_____________Resultat Förvärmning med produktflöde_______________'])
disp(['Area (m2):                                        ',num2str(A)])
disp(['Effekt (W):                                       ',num2str(q)])
disp(['Verkningsgrad:                                    ',num2str(epsilon)])
disp(['Intemperatur, kall sida (K):                      ',num2str(Tcin)])
disp(['Uttemperatur, kall sida (K):                      ',num2str(Tcut)])
disp(['Intemperatur, varm sida (K):                      ',num2str(Thin)])
disp(['Uttemperatur, varm sida (K):                      ',num2str(Thut)])
disp(['Inköpskostnad (SEK):                              ',num2str(Ka)])


%% Kylning 5-6
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan för isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan för isobuten
M_D=18.01528;    %[kg/mol]       Molmassan för vatten

%Produktflöde
%Temperaturer
Thin=711.563;                  %[K] Produktflödets temperatur in
Tut_guess=348;                 %[K] Önskad uttemp för Cp-approx, manuell iterering till denna temp nås
Tmedel=(Thin+Tut_guess)./2;    %[K] Cp tas vid medeltemp

Thut=342;

%Produktflöden [mol/s]
Fc_A=20*10^3./3600;
Fc_B=107*10^3./3600;
Fc_C=107./3600;     
Fc_D=1900*10^3./3600;

F_tot=(Fc_A+Fc_B+Fc_C+Fc_D);   %[mol/s]  Totalt molflöde

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
mvvx=20000/3600;         %[kg/s] Massflöde
Fvvx=(1000*mvvx/M_D);     %[mol/s] Molflöde
Tcin=287;                 %[K] Kondensatflödets temperatur in
cpvvx=Cp_calc(Tcin, CPcoeffDl);   %[J/(mol K)] Kondensatflödets värmekapacitet (flytande vatten)
        
U=200;                    %[W/(m2K)]      Värmegenomgångstal, gas/vätska-vvx

%Beräkning
C=[F_tot*cp_tot Fvvx*cpvvx];    %[J/(K s)]
Cmin=min(C);
Cmax=max(C);
Cmm=Cmin./Cmax;         %Kvoten Cmin/Cmax

%A=linspace(10,1000,1000);

%A=880;                      %[m2] Värmeväxlarens area

%NTU=U.*A./Cmin;

NTU=fsolve(@(NTU)NTUsolve(F_tot,cp_tot,Thut,Tcin,Cmin,Thin,NTU,Cmm),7)

epsilon=(1-exp(-NTU.*(1-Cmm)))./(1-Cmm.*exp(-NTU.*(1-Cmm)));    %Termisk verkningsgrad

%plot(A,epsilon)

A=NTU.*Cmin./U

%Kostnader och övrigt
a=32000; b=70; n=1.2;       %Kostnadsparametrar

Ka=a+b.*A.^n.*9.99;         %[SEK]  Kostnad för värmeväxlaren (1 USD = 9.99 SEK)
tdrift=8000;                %[h/år] Driftstid på ett år
Kvatten=0.05;               %[kr/kWh] Kostnad för kylvattnet

%Sökt
q=epsilon*Cmin*(Thin-Tcin);      %[W] Effekt

Kvatten_tot=(q/1000)*tdrift*Kvatten; %[SEK/år] Total kostnad för kylvattnet

Tcut=Tcin+q./(Fvvx.*cpvvx);      %[K] Uttemperatur kall sida
Thut=Thin-q./(F_tot.*cp_tot);    %[K] Uttemperatur varm sida

disp(['_____________________Resultat Kylning 5-6_____________________'])
disp(['Area (m2):                                        ',num2str(A)])
disp(['Effekt (W):                                       ',num2str(q)])
disp(['Uttemperatur, kall sida (K):                      ',num2str(Tcut)])
disp(['Uttemperatur, varm sida (K):                      ',num2str(Thut)])
disp(['Inköpskostnad (SEK):                              ',num2str(Ka)])
disp(['Kostnad för kylvattnet (SEK/år):                  ',num2str(Kvatten_tot)])


%% Ugn 1-2
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan för isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan för isobuten
M_D=18.01528;    %[kg/mol]       Molmassan för vatten

%Reaktantflöde
%Temperaturer
Thin=658.5986;           %[K]            Reaktantflödets temperatur in
Tut=929;                 %[K]            Önskad uttemperatur
Tmedel=(Thin+Tut)./2;    %[K]            Cp tas vid medeltemperatur        
dT=Tut-Thin;             %[K]            Temperaturdifferens

%Produktflöden
Fc_A=128*10^3/3600;
Fc_B=5*10^3/3600;
Fc_C=0/3600;     
Fc_D=1900*10^3/3600;          %[mol/s]

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

%Värmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);
cp_Dl=Cp_calc(Tmedel, CPcoeffDl);

%Total värmekapacitet för produktflödet
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%Produktflöde
q=cp_tot.*F_tot.*dT;     %[W] Ugnens nyttiga effekt

E=0.8;                   %Verkningsgrad ugn
qugn=(q/E)/10e6;         %[MW] Ugnens totala effekt

%Kostnader och övrigt
a=80000; b=109000; n=0.8;    %Kostnadsparametrar för cylindrisk ugn
tdrift=8000;                 %[h/år]
K=(a+b*qugn.^n)*9.99;        %[SEK]

Kel=0.30;                    %[SEK/kWh]
Kel_tot=(qugn*10e6/1000)*tdrift*Kel; %[SEK/år] Total kostnad för elen


disp(['_____________________Resultat Ugn 1-2________________________'])
disp(['Effekt (MW):                                      ',num2str(qugn)])
disp(['Intemperatur (K):                                 ',num2str(Thin)])
disp(['Uttemperatur (K):                                 ',num2str(Tut)])
disp(['Inköpskostnad (SEK):                              ',num2str(K)])
disp(['Kostnad för elen (SEK/år):                        ',num2str(Kel_tot)])

%% Ugn 3-4
clc, clear

% Givna data
% Molmassor
M_A=58.12*1000;  %[kg/mol]       Molmassan för isobutan
M_B=56.11*1000;  %[kg/mol]       Molmassan för isobuten
M_D=18.01528;    %[kg/mol]       Molmassan för vatten

%Produktflöde
%Temperaturer
Thin=905;                %[K]            Produktflödets temperatur in
Tut=950;                 %[K]            Önskad uttemperatur
Tmedel=(Thin+Tut)./2;    %[K]            Cp tas vid medeltemperatur        
dT=Tut-Thin;             %[K]            Temperaturdifferens

%Produktflöden
Fc_A=51*10^3/3600;
Fc_B=77*10^3/3600;
Fc_C=77*10^3/3600;     
Fc_D=1900*10^3/3600;  %[mol/s]

F_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/s]  Totalt molflöde

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
q=cp_tot.*F_tot.*dT;     %[W] Effekt

E=0.8;                   %Verkningsgrad ugn
qugn=(q/E)/10e6;         %[MW] Ugnens nyttiga effekt

%Kostnader och övrigt
a=80000; b=109000; n=0.8;    %Kostnadsparametrar för cylindrisk ugn
tdrift=8000;                 %[h/år]
K=(a+b*qugn.^n)*9.99;        %[SEK]

Kel=0.30;                    %[SEK/kWh]
Kel_tot=(qugn*10e6/1000)*tdrift*Kel; %[SEK/år] Total kostnad för elen


disp(['____________________Resultat Ugn 3-4________________________'])
disp(['Nyttig effekt (MW):                               ',num2str(qugn)])
disp(['Intemperatur (K):                                 ',num2str(Thin)])
disp(['Uttemperatur (K):                                 ',num2str(Tut)])
disp(['Inköpskostnad (SEK):                              ',num2str(K)])
disp(['Kostnad för elen (SEK/år):                        ',num2str(Kel_tot)])

%% Återkokare
clc,clear

Tcin_rb=304;             %[K] Produktflödets temperatur in i återkokaren
q_rb=2.813e6;            %[W] Återkokarvärme

Thin_rb=180+273;         %[K] Värmeångans temperatur in

%Molflöden [mol/s]
FcA_rb=3.33e3./3600;     %Molflöde isobutan
FcB_rb=1.67e3./3600;     %Molflöde isobuten

Fctot_rb=FcA_rb+FcB_rb;  %[mol/s] Totalt molflöde (produkt)

Fh_rb=300;               %[mol/s] Molflöde varma sidan (ånga)

%Molbråk
yA_rb=FcA_rb/Fctot_rb;
yB_rb=FcB_rb/Fctot_rb;

% Cp-koefficienter
CPcoeffD = [72.43 1.039*10^-2 -1.497*10^-6 0 ];         % VATTENÅNGA

%Värmekapacitet för ångan
Thutguess_rb=330;                         %Approx uttemp för bättre värde på Cp
Thmedel_rb=(Thutguess_rb+Thin_rb)./2;     %Ångans medeltemp
cph_rb=Cp_calc(Thmedel_rb, CPcoeffD);        %Cp för vattenångan, varma sidan

%Förångningsentalpier används istället för Cp då ingen temperaturökning
%sker på kalla sidan
dH_A=19.99e3;                             %[J/mol] Förångningsentalpi för isobutan
dH_B=20.6e3;                              %[J/mol] Förångningsentalpi för isobuten

dHtot_rb=yA_rb.*dH_A+yB_rb.*dH_B;         %[J/mol] Total ångbildningsentalpi

%Cmin, Cmax, kvot
C_rb=[Fh_rb*cph_rb Fctot_rb*dHtot_rb];
Cmin_rb=min(C_rb);
Cmax_rb=max(C_rb);
Cr_rb=Cmin_rb/Cmax_rb;

U_rb=1000;       %[W/(m2 K)] Värgmegenomgångstal, kondensor/återkokare

%Beräkning
Thut_rb=Thin_rb-q_rb./(Fh_rb.*cph_rb);
eps_rb=q_rb./(Cmin_rb.*(Thin_rb-Tcin_rb));

NTU_rb=fsolve(@(NTU)eps(eps_rb,Cr_rb,NTU),1);

A_rb=NTU_rb.*Cmin_rb./U_rb;

%Kostnader och övrigt
a_rb=32000; b_rb=70; n_rb=1.2;          %Kostnadsparametrar
K_rb=a_rb+b_rb.*A_rb.^n_rb.*9.99;       %[SEK] Inköpskostnad
tdrift=8000;                            %[h/år] Driftstid

Kanga_rb=0.16;                              %[kr/kWh] Kostnad för ångan
Kangatot_rb=(q_rb/1000)*tdrift*Kanga_rb;    %[SEK/år] Årlig kostnad för kylvattnet



disp(['_______________________Återkokare________________________'])
disp(['Area (m2):                                        ',num2str(A_rb)])
disp(['Verkningsgrad:                                    ',num2str(eps_rb)])
disp(['Effekt (W):                                       ',num2str(q_rb)])
disp(['Uttemperatur, varm sida (K):                      ',num2str(Thut_rb)])
disp(['Inköpskostnad (SEK):                              ',num2str(K_rb)])
disp(['Kostnad för ångan (SEK/år):                       ',num2str(Kangatot_rb)])


%% Kondensor
clc,clear

Thin_kd=299;             %[K] Produktflödets temperatur in i kondensorn
q_kd=2.807e6;            %[W] Kondensorvärme

Tcin_kd=273+10;          %[K] Kylvattnets temperatur in

%Molflöden varma sidan (produktflöde) [mol/s]
FhA_kd=0.175e3./3600;    %Molflöde isobutan
FhB_kd=31.7e3./3600;     %Molflöde isobuten

Fhtot_kd=FhA_kd+FhB_kd;  %[mol/s] Totalt molflöde (produkt)

Fc_kd=11000;             %[mol/s] Molflöde kalla sidan

%Molbråk
yA_kd=FhA_kd/Fhtot_kd;
yB_kd=FhB_kd/Fhtot_kd;

% Cp-koefficienter
CPcoeffDl = [32.24 0.001924 1.055e-5 -3.596e-9];   % Flytande vatten

%Värmekapacitet för kylvattnet
Tcutguess_kd=292;                         %Approx uttemp för kylvattnet
Tcmedel_kd=(Tcutguess_kd+Tcin_kd)./2;     %Kylvattnets medeltemp för Cp
cpc_kd=Cp_calc(Tcmedel_kd, CPcoeffDl);    %Cp för kylvattnet


%Förångningsentalpier används istället för Cp då ingen temperaturändring
%sker på varma sidang
dH_A=19.99e3;                             %[J/mol] Förångningsentalpi för isobutan
dH_B=20.6e3;                              %[J/mol] Förångningsentalpi för isobuten

dHtot_kd=yA_kd.*dH_A+yB_kd.*dH_B;         %[J/mol] Total ångbildningsentalpi

%Cmin, Cmax, kvot
C_kd=[Fc_kd*cpc_kd Fhtot_kd*dHtot_kd];
Cmin_kd=min(C_kd);
Cmax_kd=max(C_kd);
Cr_kd=Cmin_kd/Cmax_kd;

U_kd=1000;       %[W/(m2 K)] Värgmegenomgångstal, kondensor/återkokare

%Beräkning
Tcut_kd=Tcin_kd+q_kd./(Fc_kd.*cpc_kd);
eps_kd=q_kd./(Cmin_kd.*(Thin_kd-Tcin_kd));


NTU_kd=fsolve(@(NTU)eps(eps_kd,Cr_kd,NTU),1);

A_kd=NTU_kd.*Cmin_kd./U_kd;

%Kostnader och övrigt
a_kd=32000; b_kd=70; n_kd=1.2;              %Kostnadsparametrar
K_kd=a_kd+b_kd.*A_kd.^n_kd.*9.99;           %[SEK] Inköpskostnad
tdrift=8000;                                %[h/år] Driftstid

Kvatten_kd=0.05;                                   %[SEK/kWh] Kostnad för kylvattnet (under 14 C)
Kvattentot_kd=(q_kd/1000)*tdrift*Kvatten_kd;    %[SEK/år] Total kostnad för kylvattnet



disp(['_______________________Kondensor_________________________'])
disp(['Area (m2):                                        ',num2str(A_kd)])
disp(['Verkningsgrad:                                    ',num2str(eps_kd)])
disp(['Effekt (W):                                       ',num2str(q_kd)])
disp(['Uttemperatur, kall sida (K):                      ',num2str(Tcut_kd)])
disp(['Inköpskostnad (SEK):                              ',num2str(K_kd)])
disp(['Kostnad för kylvattnet (SEK/år):                  ',num2str(Kvattentot_kd)])


%% Funktion för beräkning av Cp(T)
function Cp = Cp_calc(T, CPcoeff)
    A = CPcoeff(1);
    B = CPcoeff(2);
    C = CPcoeff(3);
    D = CPcoeff(4);

    Cp = A + B*T + C*T.^2 + D*T.^3;

end

%% Beräkning av NTU eller epsilon
function diff=eps(epsilon,Cmm,NTU)

VL=epsilon;
HL=(1-exp(-NTU.*(1-Cmm)))./(1-Cmm.*exp(-NTU.*(1-Cmm)));

diff=VL-HL;

end

%% AAAAA
function diff2=NTUsolve(F_tot,cp_tot,Thut,Tcin,Cmin,Thin,NTU,Cmm)

e1=-F_tot.*cp_tot.*(Thut-Thin)./(Cmin.*(Thin-Tcin));

e2=(1-exp(-NTU.*(1-Cmm)))./(1-Cmm.*exp(-NTU.*(1-Cmm)));

diff2=e1-e2;

end