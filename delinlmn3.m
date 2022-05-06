% Reaktor KAA146 Grundläggande kemiteknik, Projektarbete grupp 4
clc, clear, close all

% isobutan = A
% isobuten = B
% vätgas   = C
% vatten   = D

% DATA
% Molmassor
MA = 58.1222;      % g mol-1  Butan
MB = 56.1063;      % g mol-1  Buten
MC = 2.0160;       % g mol-1  Vätgas
MD = 18.0160;      % g mol-1  Vatten

R = 8.31447;        % Gaskonstanten [J mol^-1 K^-1]
P = 1;              % [bar]
T_reaktor1 = 950;   % [K]
T_reaktor2 = 950;   % [K]
rho_cat = 1120;     % [kg m3^-1]
Ea = 141e3;         % [J mol^-1]
k = 0.0596;         % [mol kg cat.^-1 s^-1 bar^-1 vid 550 C]
A = k/(exp(-Ea/(R*(550+273.15))));
K1 = 22.9;          % [bar^-1]
K2 = 7.56;          % [bar^-1]

CPcoeffA = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];     % ISOBUTAN
CPcoeffB = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];     % ISOBUTEN
CPcoeffC = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];    % H2
CPcoeffD = [32.24 0.001924 1.055*10^-5 -3.596*10^-9];   % H2O
CPcoeffDl = [72.43 1.039*10^-2 -1.497*10^-6 0 ];         % H2O flytande

Cp = [CPcoeffA; CPcoeffB; CPcoeffC; CPcoeffD; CPcoeffDl];


%% REAKTOR
dH0A   = -134.2e3;       % [kJ mol^-1]
dH0B   = -17.9e3;        % [kJ mol^-1]
dH0C   = 0;              % [kJ mol^-1]
dHr0   = dH0B+dH0C-dH0A; % [kJ mol^-1]

% REAKTOR 1
U01 = [133/3.6 0/3.6 0 1330/3.6 T_reaktor1]; %[mol s^-1]

Wstart1 = 0; %Massa cat. [kg]
Wfinal1 = 3000; %[kg]
Wspan1 = [Wstart1 Wfinal1];
[W1,U1] = ode15s(@PBR_ode,Wspan1,U01,[],dHr0,A,Ea,K1,K2,P,Cp); % Löser molflöden för reaktor 1

T1 = U1(:,5); %[K]
% V1= W1./rho_cat; %[m3]
X1 = (U1(1,1)-U1(end,1))./U1(1,1);

% REAKTOR 2
U02 = [U1(end,1) U1(end,2) U1(end,3) U1(end,4) T_reaktor2];
Wstart2 = 0; %Massa cat. [kg]
Wfinal2 = 4000; 
Wspan2 = [Wstart2 Wfinal2];
[W2,U2] = ode15s(@PBR_ode,Wspan2,U02,[],dHr0,A,Ea,K1,K2,P,Cp); % Löser molflöden för reaktor 2

T2 = U2(:,5); %[K]
% V2 = W2./rho_cat; %[m3]
X2 = (U2(1,1)-U2(end,1))./U2(1,1);

% Alternativ använda en reaktor
Wstart_alt1 = 0; %Massa cat. [kg]
Wfinal_alt1 = Wfinal1+Wfinal2; 
Wspan_alt1 = [Wstart_alt1 Wfinal_alt1];
[W_alt1,U_alt1] = ode15s(@PBR_ode,Wspan_alt1,U01,[],dHr0,A,Ea,K1,K2,P,Cp);

T_alt = U_alt1(:,5); %[K]
V_alt= W_alt1./rho_cat; %[m3]

X_alt1= (U_alt1(1,1)-U_alt1(:,1))/U_alt1(1,1);

% BÅDA REAKTORERNA
U = [U1; U2(2:end,:)];
X = (U(1,1)-U(:,1))/U(1,1);
W = [W1; W2(2:end,:)+W1(end)]; % [kg]
T = [T1; T2(2:end,:)];           % [K]

% Molflöden ut ur reaktor 2
FA = U(:,1);
FB = U(:,2);
FC = U(:,3); 


% PLOTTAR:
% Plottar omsättningsgraden mot massan katalysator för de båda reaktorerna
figure(1);
plot(W_alt1,X_alt1,'r--')
hold on
plot(W,X,'Color','#0072BD')
title('Reaktor 1 och 2'),xlabel('Katalysatormassa, W [kg]'), ylabel('Omsättningsgrad, X')
hold on
plot([Wfinal1 Wfinal1],[0 1],'k--')     % Avgränsning mellan reaktor 1 och 2

% Plottar temperaturen mot omsättningsgraden för de båda reaktorerna
figure(2);
plot(X,T)
title('Reaktor 1 och 2'),xlabel('Omsättningsgrad, X'),ylabel('T [K]')

% Plottar molflöden för de olika specierna mot katalysatormassan
figure(3);
plot(W,FA,W,FB)
hold on
plot([Wfinal1 Wfinal1],[0 40],'k--')   % Avgränsning mellan reaktor 1 och 2
title('Molflöde för specier vs katalysatormassa')
xlabel('Katalysatormassa, W [kg]'), ylabel('Molflöde [mol/s]')
legend('isobutan','isobuten och vätgas','Location','northwest')


%% Tryckkärlsväggens (reaktorväggens) tjocklek
% REKTOR 1
T_F1 = (9/5)*(T_reaktor1-273.15)+32;	                  % Temperatur i Farenheit, om över 900F måste rostfritt stål användas, max 1500F
Vol1 = Wfinal1./rho_cat;
D1 = (2*Vol1/pi)^(1/3);                                   % Diameter på reaktor 1 [m]
Smax = 74.5E6;                                            % Maximalt tillåtna spänningen över 900F [N/m2]
P_konstr = P*10^5*1.1;                                    % Konstruktionstryck, 10% mer än arbetstryck [Pa]
E = 1;                                                    % Svetsverkningsgrad
wall1= ((P_konstr*D1)/((2*Smax*E)-(1.2*P_konstr))).*10^3; %[mm]
rho_wall = [7900 8000];

if D1 < 1
    wall1 = 6;
elseif D1 >= 1 && D1 < 2
        wall1 = 7;
elseif D1 >= 2 && D1 < 2.5
     wall1 = 9;
end

V_wall1 = pi.*((D1+2.*wall1*10^-3)/2).^2.*2*(D1+2*wall1*10^-3) - pi*(D1/2)^2*2*D1; %[m3]
mass_wall1 = V_wall1(1)*rho_wall(1);    %[kg]

% Horisontell reaktor pris
Reaktor_param = [12800 73 0.85];

cost_reak1 = Cost(mass_wall1,Reaktor_param);        % Kostnad för reaktor 1 år 2010 (CEPCI = 532.9)
cost_reak1_cat = cost_reak1*1.5;                    % Kostnad för reaktor 1 + katalysator år 2010 (CEPCI = 532.9)
cost1_2020 = (cost_reak1_cat*(569/532.9))*9.99*4;   % Kostnad för reaktor 1 + katalysator samt montering år 2020 i SEK

% REKTOR 1 alt.
T_Falt = (9/5)*(T_reaktor1-273.15)+32;                           % Temperatur i Farenheit, om över 900F måste rostfritt stål användas, max 1500F
Vol_alt = Wfinal_alt1./rho_cat;                                  % Volym på reaktor 1 (alternativ) [m^3]
D_alt = (2*Vol_alt/pi)^(1/3);                                    % Diameter på reaktor 1 (alternativ) [m]
wall_alt = ((P_konstr*D_alt)./((2*Smax*1)-(1.2*P_konstr))).*10^3; %[mm]

if D_alt < 1
    wall_alt = 6;
elseif D_alt >= 1 && D_alt < 2
   wall_alt = 7;
elseif D_alt >= 2 && D_alt < 2.5
   wall_alt = 9;
end

V_wall_alt = pi.*((D_alt+2.*wall_alt*10^-3)/2).^2.*2*(D_alt+2*wall_alt*10^-3) - pi*(D_alt/2)^2*2*D_alt; %[m3]
mass_wall_alt = V_wall_alt.*rho_wall(1);    %[kg]

% Horisontell reaktor pris
cost_reak_alt = Cost(mass_wall_alt,Reaktor_param);
cost_reak_alt_cat = cost_reak_alt*1.5;                       % Kostnad för reaktor 1 + katalysator år 2010 (CEPCI = 532.9)
cost_alt_2020 = (cost_reak_alt_cat*(569/532.9))*9.99*4;      % Kostnad för reaktor 1 + katalysator samt montering år 2020 i SEK

% REAKTOR 2
T_F2 = (9/5)*(T_reaktor2-273.15)+32;                        % Temperatur i Farenheit, om över 900F måste rostfritt stål användas, max 1500F
Vol2 = Wfinal2./rho_cat;
D2 = (2*Vol2/pi)^(1/3);                                     % Diameter på reaktor  [m]
wall2= ((P_konstr*D2)./((2*Smax*1)-(1.2*P_konstr))).*10^3;  % [mm]

if D2 < 1
    wall2 = 6;
elseif D2 >= 1 && D2 < 2
    wall2 = 7;
elseif D2 >= 2 && D2 < 2.5
    wall2 = 9;
end

V_wall2 = pi.*((D2+2.*wall2*10^-3)/2).^2.*2*(D2+2*wall2*10^-3) - pi*(D2/2)^2*2*D2; %[m3]
mass_wall2 = V_wall2.*rho_wall(1);  %[kg]

% Horisontell reaktor
cost_reak2 = Cost(mass_wall2,Reaktor_param);
cost_reak2_cat = cost_reak2*1.5;                         % Kostnad för reaktor 2 + katalysator år 2010 (CEPCI = 532.9)
cost2_2020 = (cost_reak2_cat*(569/532.9))*9.99*4;        % Kostnad för reaktor 2 + katalysator samt montering år 2020 i SEK
cost_allareakt = cost1_2020 + cost2_2020;

%% Utskrivning av resultat: Reaktor
disp('REAKTOR:')
disp(' ' )

% Ut ur reaktor 1
disp('______________________Reaktor 1__________________________')
disp(['Butan ut (mol/s):                                ',num2str(U1(end,1))])
disp(['Buten ut (mol/s):                                ',num2str(U1(end,2))])
disp(['Vätgas ut (mol/s):                               ',num2str(U1(end,3))])
disp(['Temperatur ut (K):                               ',num2str(U1(end,5))])
disp(['Katalysatormassa (kg):                           ',num2str(Wfinal1)])
disp(['Volym (m^3):                                     ',num2str(Vol1)])
disp(['Diameter (m):                                    ',num2str(D1)])
disp(['Längd (m):                                       ',num2str(D1*2)])
disp(['Omsättningsgrad:                                 ',num2str(X1)])
disp(' ')

% Ut ur reaktor 2
disp('______________________Reaktor 2__________________________')
disp(['Butan ut (mol/s):                                ',num2str(U(end,1))])
disp(['Buten ut (mol/s):                                ',num2str(U(end,2))])
disp(['Vätgas ut (mol/s):                               ',num2str(U(end,3))])
disp(['Temperatur ut (K):                               ',num2str(U(end,5))])
disp(['Katalysatormassa (kg):                           ',num2str(Wfinal2)])
disp(['Volym (m^3):                                     ',num2str(Vol2)])
disp(['Diameter (m):                                    ',num2str(D2)])
disp(['Längd (m):                                       ',num2str(D2*2)])
disp(['Omsättningsgrad:                                 ',num2str(X2)])

disp(' ')
disp(['Omsättningsgrad TOTAL:                           ',num2str(X(end))])
disp(['Katalysatormassa (kg) TOTAL:                     ',num2str(Wfinal1+Wfinal2)])
disp(['Volym (m^3) TOTAL:                               ',num2str(Vol1+Vol2)])
disp(' ')

% Reaktorkostnad
disp('______________________Reaktorkostnader________________________')
disp(['Kostnad för reaktor 1 år 2020 (SEK):              ',num2str(round(cost1_2020, -3))])
disp(['Kostnad för alternativ reaktor 1 år 2020 (SEK):   ',num2str(round(cost_alt_2020, -3))])
disp(['Kostnad för reaktor 2 år 2020 (SEK):              ',num2str(round(cost2_2020, -3))])
disp(['Kostnad för reaktor 1 och 2 år 2020 (SEK):        ',num2str(round(cost_allareakt, -3))])
disp('______________________________________________________________')
disp(' ')

%% Värmeväxling
%förvärmning
%Temperaturer
Thin=U(end,5);        %[K]        Produktflödets intemperatur (varm sida)
Tcin=180+273;    %[K]        Reaktantflödets intemperatur (kall sida) (blandat med 180C ånga)
Tmedel=(Thin+Tcin)./2; %[K]  Cp tas vid medeltemp

%Flöde varm sida [mol/s]
Fh_A=U(end,1);
Fh_B=U(end,2);
Fh_C=U(end,3);     
Fh_D=U(end,4);

Fh_tot=(Fh_A+Fh_B+Fh_C+Fh_D);  %[mol/s]  Totalt molflöde varm sida

%Flöden kall sida
Fc_A=U01(1);
Fc_B=U01(2);
Fc_C=U01(3);     
Fc_D=U01(4);

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

%Värmekapaciteter [J/(mol K)]
cp_A=Cp_calc(Tmedel, CPcoeffA);        
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);

%Total värmekapacitet, varm sida
cph_tot=yh_A.*cp_A+yh_B.*cp_B+yh_C.*cp_C+yh_D.*cp_D; %[J/(molK)] 

%Total värmekapacitet, kalla sida
cpc_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

U_VVX1=50;                   %[W/(m2 K)]      Värmegenomgångstal, gas/gas-vvx
A_VVX1=1000;                  %[m2] Area
%A=linspace(0,1000,1000)    %För att plotta epsilon mot A

%Kostnader
a=32000; b=70; n=1.2;   %Kostnadsparametrar
Ka=a+b.*A_VVX1.^n.*9.99;     %[SEK]    Kostnad för värmeväxlaren (1 USD = 9.99 SEK)

%Beräkning
C=[Fh_tot*cph_tot Fc_tot*cpc_tot];
Cmin=min(C);
Cmax=max(C);
Cmm=Cmin./Cmax;         %Kvoten Cmin/Cmax

NTU=U_VVX1*A_VVX1/Cmin;           %Number of transfered units


epsilon=(1-exp(-NTU.*(1-Cmm)))./(1-Cmm.*exp(-NTU.*(1-Cmm)));    %Termisk verkningsgrad


%Sökt
q=epsilon*Cmin*(Thin-Tcin);          %[W] Effekt
Tcut=Tcin+q./(Fc_tot.*cpc_tot);      %[K] Uttemperatur kall sida
Thut=Thin-q./(Fh_tot.*cph_tot);      %[K] Uttemperatur varm sida

%plot(A,epsilon) %För optimering av area

disp('_____________Resultat Förvärmning med produktflöde_______________')
disp(['Area (m2):                                        ',num2str(A_VVX1)])
disp(['Effekt (MW):                                      ',num2str(q*1e-6)])
disp(['Verkningsgrad:                                    ',num2str(epsilon)])
disp(['Intemperatur, kall sida (K):                      ',num2str(Tcin)])
disp(['Uttemperatur, kall sida (K):                      ',num2str(Tcut)])
disp(['Intemperatur, varm sida (K):                      ',num2str(Thin)])
disp(['Uttemperatur, varm sida (K):                      ',num2str(Thut)])
disp(['Kostnad för värmeväxlaren (SEK):                  ',num2str(Ka)])

%Ugn 1-2
%Temperaturer
Tcin=Tcut;               %[K]            Reaktantflödets temperatur in
Tut=950;                 %[K]            Önskad uttemperatur
Tmedel=(Tcin+Tut)./2;    %[K]            Cp tas vid medeltemperatur        
dT=Tut-Tcin;             %[K]            Temperaturdifferens

%Värmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);

%Total värmekapacitet för produktflödet
cp_tot=yc_A.*cp_A+yc_B.*cp_B+yc_C.*cp_C+yc_D.*cp_D; %[J/(molK)] 

%Produktflöde
q=cp_tot.*Fc_tot.*dT;     %[W] Ugnens nyttiga effekt

E=0.8;                   %Verkningsgrad ugn
qugn=(q/E)/10e6;         %[MW] Ugnens totala effekt

%Kostnader och övrigt
a=80000; b=109000; n=0.8;    %Kostnadsparametrar för cylindrisk ugn
tdrift=8000;                 %[h/år]
K=(a+b*qugn.^n)*9.99;        %[SEK]

Kel=0.30;                    %[SEK/kWh]
Kel_tot=(qugn*10e6/1000)*tdrift*Kel; %[SEK/år] Total kostnad för elen


disp('________________________Resultat Ugn 1-2___________________________')
disp(['Effekt (MW):                                      ',num2str(qugn)])
disp(['Intemperatur (K):                                 ',num2str(Tcin)])
disp(['Uttemperatur (K):                                 ',num2str(Tut)])
disp(['Kostnad för ugnen (SEK):                          ',num2str(K)])
disp(['Kostnad för elen (SEK/år):                        ',num2str(Kel_tot)])

%Ugn 3-4
Tcin=U1(end, 5);         %[K]            Produktflödets temperatur in
Tut=950;                 %[K]            Önskad uttemperatur
Tmedel=(Tcin+Tut)./2;    %[K]            Cp tas vid medeltemperatur        
dT=Tut-Tcin;             %[K]            Temperaturdifferens

%Produktflöden
Fc_A= U1(end, 1);
Fc_B= U1(end, 2);
Fc_C= U1(end, 3);     
Fc_D= U1(end, 4);  %[mol/s]

F_tot=Fc_A+Fc_B+Fc_C+Fc_D;  %[mol/s]  Totalt molflöde

%Molbråk
yc_A=Fc_A./F_tot;
yc_B=Fc_B./F_tot;
yc_C=Fc_C./F_tot;
yc_D=Fc_D./F_tot;

%Värmekapaciteter
cp_A=Cp_calc(Tmedel, CPcoeffA);        %[J/(molK)]
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);

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


disp('____________________Resultat Ugn 3-4________________________')
disp(['Nyttig effekt (MW):                               ',num2str(qugn)])
disp(['Intemperatur (K):                                 ',num2str(Tcin)])
disp(['Uttemperatur (K):                                 ',num2str(Tut)])
disp(['Kostnad för ugnen (SEK):                          ',num2str(K)])
disp(['Kostnad för elen (SEK/år):                        ',num2str(Kel_tot)])


%Kylning 6-7
%Produktflöde
%Temperaturer
Thin=U(end, 5);                %[K] Produktflödets temperatur in
Tut_guess=348;                 %[K] Önskad uttemp för Cp-approx, manuell iterering till denna temp nås
Tmedel=(Thin+Tut_guess)./2;    %[K] Cp tas vid medeltemp

%Produktflöden [mol/s] (varm sida)
Fh_A=20*10^3./3600;
Fh_B=107*10^3./3600;
Fh_C=107./3600;     
Fh_D=1900*10^3./3600;

Fh_tot=(Fh_A+Fh_B+Fh_C+Fh_D);   %[mol/s]  Totalt molflöde

%Molbråk
yh_A=Fh_A./Fh_tot;
yh_B=Fh_B./Fh_tot;
yh_C=Fh_C./Fh_tot;
yh_D=Fh_D./Fh_tot;

%Värmekapaciteter [J/(mol K)]
cp_A=Cp_calc(Tmedel, CPcoeffA);        
cp_B=Cp_calc(Tmedel, CPcoeffB);
cp_C=Cp_calc(Tmedel, CPcoeffC);
cp_D=Cp_calc(Tmedel, CPcoeffD);

%Total värmekapacitet för flödet
cp_tot = yh_A.*cp_A + yh_B.*cp_B + yh_C.*cp_C + yh_D.*cp_D; %[J/(molK)] 

%VVX-sidan (kalla)
mvvx=23350/3600;                  %[kg/s] Massflöde
Fvvx=(1000*mvvx/MD);              %[mol/s] Molflöde
Tcin=287;                         %[K] Kondensatflödets temperatur in
cpvvx=Cp_calc(Tcin, CPcoeffDl);   %[J/(mol K)] Kondensatflödets värmekapacitet (flytande vatten)
        
U_VVX2=200;                            %[W/(m2K)] Värmegenomgångstal, gas/vätska-vvx

%Beräkning
C=[Fh_tot*cp_tot Fvvx*cpvvx];    %[J/(K s)]
Cmin=min(C);
Cmax=max(C);
Cmm=Cmin./Cmax;                  %Kvoten Cmin/Cmax

%A=linspace(10,1000,1000);

A_VVX2=880;                           %[m2] Värmeväxlarens area

NTU=U_VVX2.*A_VVX2./Cmin;

epsilon=(1-exp(-NTU.*(1-Cmm)))./(1-Cmm.*exp(-NTU.*(1-Cmm)));    %Termisk verkningsgrad

%plot(A,epsilon)

%Kostnader och övrigt
a=32000; b=70; n=1.2;       %Kostnadsparametrar

Ka=a+b.*A_VVX2.^n.*9.99;         %[SEK]  Kostnad för värmeväxlaren (1 USD = 9.99 SEK)
tdrift=8000;                %[h/år] Driftstid på ett år
Kvatten=0.05;               %[kr/kWh] Kostnad för kylvattnet

%Sökt
q=epsilon*Cmin*(Thin-Tcin);      %[W] Effekt

Kvatten_tot=(q/1000)*tdrift*Kvatten; %[SEK/år] Total kostnad för kylvattnet

Tcut=Tcin+q./(Fvvx.*cpvvx);      %[K] Uttemperatur kall sida
Thut=Thin-q./(Fh_tot.*cp_tot);   %[K] Uttemperatur varm sida

disp(' ')
disp('_____________________Resultat Kylning 6-7_____________________')
disp(['Area (m2):                                        ',num2str(A_VVX2)])
disp(['Effekt (MW):                                      ',num2str(q*1e-6)])
disp(['Verkningsgrad:                                    ',num2str(epsilon)])
disp(['Intemperatur, kall sida (K):                      ',num2str(Tcin)])
disp(['Uttemperatur, kall sida (K):                      ',num2str(Tcut)])
disp(['Intemperatur, varm sida (K):                      ',num2str(Thin)])
disp(['Uttemperatur, varm sida (K):                      ',num2str(Thut)])
disp(['Kostnad för värmeväxlaren (SEK):                  ',num2str(Ka)])
disp(['Kostnad för kylvattnet (SEK/år):                  ',num2str(Kvatten_tot)])


%% Flash
% Data
k = 0.107; % m s-1 @8bar
P = 101325; % Pa (1atm)
T = 348; % K

%Antoinekonstanter A  B  C
Ant1 =  [15.7564 2132.42 -33.15];  % buten
Ant2 =  [15.6782 2154.90 -34.42];  % butan 
%Wilsonfaktorer
W12 = 0.48584; 
W21 = 1.64637; 

rhoA = 559.0;      % kg m-3 Buten
rhoB = 556.62;     % kg m-3 Butan
rhoD = 974.93;     % kg m-3 Vatten @ 75C

% Vätskeflöde från tanken (Vatten)
L1 = (U(end,4)*MD*1e-3)/rhoD; %m3/s
% Ångflöde ut
V1 = ((sum(U(end,:))-U(end,4)-U(end,5))*R*T)/P; %m3/s

%uppehållstid
tau = 10*60; % s (10min)
rho_L1 = rhoD;
rho_V1 = P*1e-3*(U(end,1)*MA + U(end,2)*MB + U(end,3)*MC)/(T*R*V1); %kg m-3
ut = k*sqrt((rho_L1 - rho_V1)/rho_V1); %m s-1

% Diameter
D1 = sqrt(4*V1/(pi*0.15*ut)); %m

% Vätskehöjd
HL1 = (L1*tau)/((pi*D1.^2)/4); %m

%Höjd flashtank
H1 = HL1 + 1.5*D1; %m

% Kärlets tjocklek tank 1
rho_wall = [7900 8000];            % Densitet ( kolstål / rostfritt stål )
S = [88.9 120.65]*10^6;            % Materialspänning [N m-1]

t_flash1 = (1.1*P*D1*10^3)./(2*S-1.2*1.1*P); % [mm]

if t_flash1(1) < 12 && t_flash1(2) < 12
    t_flash1 = [12, 12];
end

Vwall_flash1 = pi.*((D1+2.*t_flash1*10^-3)/2).^2.*(H1+2*t_flash1*10^-3) - pi*(D1/2)^2*H1;
mwall_flash1 = Vwall_flash1.*rho_wall;


%% tank 2
P = 101325*6; % Pa (6atm)

%molbråk
x1 = U(end,2)/(U(end,1) + U(end,2));

%BP-beräkning
options = optimset('Display','off');    % Så att skit inte skrivs ut efter fsolve
[gamma1, gamma2] = wilson(x1,W12,W21);
Tguess = 273;
T = fsolve(@(T) find_Tb(T,x1,gamma1,gamma2,Ant1,Ant2,P*0.0075006168), Tguess, options);

%Ång- och vätskeflöde ut ur tanken
V2 = (U(end,3)*R*T)/P; %m3/s
L2 = (U(end,1)*MA*1e-3)/rhoA + (U(end,2)*MB*1e-3)/rhoB; %m3/s

%Densiteter för ånga och vätska
rho_L2 = (U(end,1)*rhoA + U(end,2)*rhoB)/(U(end,1) + U(end,2)); %kg m-3
rho_V2 = P*1e-3*(U(end,3)*MC)/(T*R*V2); %kg m-3

ut = k*sqrt((rho_L2 - rho_V2)/rho_V2); %m s-1

% Diameter
D2 = sqrt(4*V2/(pi*0.15*ut)); %m

% Vätskehöjd
HL2 = (L2*tau)/((pi*D2.^2)/4); %m

%Höjd flashtank
H2 = HL2 + 1.5*D2; %m

% Kärlets tjocklek tank 2
S = [88.9 137.9]*10^6;    % Materialspänning [N m-1]

t_flash2 = (1.1*P*D1*10^3)./(2*S-1.2*1.1*P); % [mm]

Vwall_flash2 = pi.*((D2+2.*t_flash2*10^-3)/2).^2.*(H2+2*t_flash2*10^-3) - pi*(D2/2)^2*H2;
mwall_flash2 = Vwall_flash2.*rho_wall;

%% Kostnader flash
kurs = 9.99;                        % Växelkursen sek/dollar
lang = 4;                           % Langfaktorn
index = 596/532.9;

Param_skalmassa = [11600 34 0.85
                   17400 79 0.85];
             
% Flash 1               
kostnad_flash1_kol = Cost(mwall_flash1(1),Param_skalmassa(1,:))*kurs*lang*index;
kostnad_flash1_rostfri = Cost(mwall_flash1(2),Param_skalmassa(2,:))*kurs*lang*index;

% Flash 2              
kostnad_flash2_kol = Cost(mwall_flash2(1),Param_skalmassa(1,:))*kurs*lang*index;
kostnad_flash2_rostfri = Cost(mwall_flash2(2),Param_skalmassa(2,:))*kurs*lang*index;

%% UTSKRIVNING AV RESULTAT: Separation (flashtankar)
disp(' ')
disp('SEPARATION (flashtankar):')
disp(' ')
disp('______________________Dimensionering__________________________')
disp(['Diameter tank 1 (m):                 ' num2str(D1, '%.2f')])
disp(['Höjd på tank 1 (m):                  ' num2str(H1, '%.2f')])
disp(' ')
disp(['Diameter tank 2 (m):                 ' num2str(D2, '%.2f')])
disp(['Höjd på tank 2 (m):                  ' num2str(H2, '%.2f')])
disp(' ')
disp('______________________Utrustningskostnader_____________________')
disp('Tank 1:')
disp(['Väggtjocklek (kolstål) (mm):        ' num2str(t_flash1(1), '%.2f')])
disp(['Väggtjocklek (rostfritt stål) (mm): ' num2str(t_flash1(2), '%.2f')])
disp(['Volym (kolstål) (m^3)               ' num2str(Vwall_flash1(1), '%.2f')])
disp(['Volym (rostfritt stål) (m^3)        ' num2str(Vwall_flash1(2), '%.2f')])
disp(['Massa (kolstål) (kg)                ' num2str(mwall_flash1(1), '%.2f')])
disp(['Massa (rostfritt stål) (kg)         ' num2str(mwall_flash1(2), '%.2f')])
disp(['Kostnad (kolstål) (kr)              ' num2str(kostnad_flash1_kol, '%.2f')])
disp(['Kostnad (rostfritt stål) (kr)       ' num2str(kostnad_flash1_rostfri, '%.2f')])
disp(' ')
disp('Tank 2:')
disp(['Temperatur:                         ' num2str(T)])
disp(['Väggtjocklek (kolstål) (mm):        ' num2str(t_flash2(1), '%.2f')])
disp(['Väggtjocklek (rostfritt stål) (mm): ' num2str(t_flash2(2), '%.2f')])
disp(['Volym (kolstål) (m^3)               ' num2str(Vwall_flash2(1), '%.2f')])
disp(['Volym (rostfritt stål) (m^3)        ' num2str(Vwall_flash2(2), '%.2f')])
disp(['Massa (kolstål) (kg)                ' num2str(mwall_flash2(1), '%.2f')])
disp(['Massa (rostfritt stål) (kg)         ' num2str(mwall_flash2(2), '%.2f')])
disp(['Kostnad (kolstål) (kr)              ' num2str(kostnad_flash2_kol, '%.2f')])
disp(['Kostnad (rostfritt stål) (kr)       ' num2str(kostnad_flash2_rostfri, '%.2f')])

%% Destillation

% Data - separation 
q = 1;              % kokvarmt tillflöde
F = (U(end,1) + U(end,2))*3.6;     % kmol h-1
P = 2280;           % mmHg (3 atm)
xf = U(end,2)/(U(end,1)+U(end,2)); % Komposition in till destillationskolonnen, molbråk buten
xd = 0.95;          % destillatbråk 
xb = 0.05;          % bottenbråk
Rec = 3;           % återflödesförhållande

% Kokpunkter @ 3 atm
Tb_1 = -6.3+273.15; % K   (buten)
Tb_2 = -1+273.15;   % K   (butan)

%Densiteter
L_rho1 = 559.0;          % kg m-3 Buten
L_rho2 = 556.62;         % kg m-3 Butan
%Ytspänning
surfaceten = 24;          % dyn cm-1
    
% Beräkning av flöden 
B = F*((xf-xd)/(xb-xd));
D = F-B;

L = Rec*D; 
V = L+D;

l = L + q*F;    % L-streck
v = l-B;        % V-streck

% Återkokare 
T = zeros(1,30);
[gamma1, gamma2] = wilson(xb,W12,W21);
Tstart = 273.15;
TB=fsolve(@(T)find_Tb(T,xb,gamma1,gamma2,Ant1,Ant2,P),Tstart,options);
TB_reboiler = TB;
Psat1 = antoine(TB, Ant1);
y0 = (gamma1*xb.*Psat1)/P;
x1 = (v/l)*y0 + (L_rho1/l)*xb;

% Avdrivardel 
x = zeros(1,60);
y = zeros(1,60);
x(1) = x1;
y(1)=y0;
i = 0; 

% Sorels metod
% Avdrivare
while x<xf
    i = i+1; 
    [gamma1, gamma2] = wilson(x(i),W12,W21);
    Tstart = TB;
    TB=fsolve(@(T)find_Tb(T,x(i),gamma1,gamma2,Ant1,Ant2,P),Tstart,options);
    T(i) = TB;
    Psat1 = antoine(TB, Ant1);                                         % Ångtryck
    y(i) = (gamma1*x(i).*Psat1)/P;
    x(i+1) = (v/l)*y(i) + (B/l)*xb;
end

m = i+1;

% Förstärkare
while y(i)<xd
    x(i+1)= (V/L)*y(i) + (1/L)*(B*xb-F*xf);                        % Komponentbalans över förstärkardelen
    i=i+1;
    [gamma1, gamma2] = wilson(x(i),W12,W21);
    Tstart = TB;
    TB=fsolve(@(T)find_Tb(T,x(i),gamma1,gamma2,Ant1,Ant2,P),Tstart,options);
    T(i) = TB;
    Psat1 = antoine(TB, Ant1);                                      % Ångtryck
    y(i) = (gamma1*x(i).*Psat1)/P;
end

bottnar_ideal = i;
bottnar_verklig = round(i/0.7);

% Jämviktkurva och jämviktsplot
xeq = 0:0.001:1;    
[gamma1, gamma2] = wilson(xeq,W12,W21);
Tstart = linspace(-6.3+273.15,-1+273.15,1001);
TBeq=fsolve(@(T)find_Tb(T,xeq,gamma1,gamma2,Ant1,Ant2,P),Tstart,options);
Psat1 = antoine(TBeq, Ant1);                                            % Ångtryck
yeq = (gamma1.*xeq.*Psat1)./P;

figure(4);
plot(xeq,xeq)   % Referenslinje
hold on
plot(xeq,yeq)
xlabel('x_1'), ylabel('y_1')
axis([0 1 0 1])

legend('Referenslinje', 'Jämviktskurva','Location','northwest')

%% Dimensionering

% Sammansättning ut ur återkokare
v_x1 = y0;              % Buten
v_x2 = 1-y0;            % Butan
% Sammansättning in i återkokare
l_x1 = x(1);
l_x2 = 1 - x(1);

% Flödesfaktorer
surftention = 15.28;        %dyn cm-1
Fst = (surftention/20)^0.2;
vaporveloc = 0.7;
Ff = 1; %nono-foaming
Fha = 1; %Hålen är bra

%Bottenavstånd
trayheight = 0.45; %m

%Densiteter för vätska och gas
rho_L = ((l_x1*MB*1e-3)/(l_x1*MB*1e-3 + l_x2*MA*1e-3))*L_rho1 + ((l_x2*MA*1e-3)/(l_x1*MB*1e-3 + l_x2*MA*1e-3))*L_rho2; % kg m-3
rho_V = v_x1*MB*1e-3*(P*133.322368/(R*TB_reboiler)) + v_x2*MA*1e-3*(P*133.322368/(R*TB_reboiler));

%molmassor
M_L = l_x1*MB + l_x2*MA;
M_V = v_x1*MB + v_x2*MA;

%Belastningsparameter
Flv = ((l*M_L)/(v*M_V)) * sqrt(rho_V/rho_L);

% från diagramet
Cf = 0.25*0.3048; %ft/s -> m/s

%Flödningsparametern
C = Fst*Ff*Fha*Cf;

% Ånghastigheten vid flödning
Uf = C*sqrt((rho_L-rho_V)/rho_V);

%Ånghastighet
Uvap = vaporveloc*Uf;

%Aktiv area
Aaktiv = (V*(1/3.6)*M_V*1e-3*(1/rho_V))/Uvap;

%Total area
Atot = Aaktiv/0.8;

%Bottendiameter
d = 2*sqrt(Atot/pi);

%Kolonnens höjd
h = trayheight * (bottnar_verklig + 1);

% Kärlets tjocklek
rho_wall = [7900 8000];            % Densitet ( kolstål / rostfritt stål )
t = 9;  % [mm]

Vwall_dest = pi.*((d+2.*t*10^-3)/2).^2.*(h+2*t*10^-3) - pi*(d/2)^2*h;
mwall_dest = Vwall_dest.*rho_wall;

% Värmen
Hvap1 = 20.6e3;  % j mol-1
Hvap2 = 19.99e3; % j mol-1

% Värmen
Q_condensor = (Hvap1*V*1e3*x(end) + Hvap2*V*1e3*(1-x(end))) * 3600^-1;    % W
Q_reboiler = (Hvap1*v*1e3*y0 + Hvap2*v*1e3*(1-y0)) * 3600^-1;             % W

%% Kostnader: Separation

% Bottnar
% Parametrar givna i PM
Param_bottnar = [130 440 1.8        % sieve tray
                 210 400 1.9        % valve tray
                 340 400 1.9];      % bubble cap tray
     
kurs = 9.99;                        % Växelkursen sek/dollar
lang = 4;                           % Langfaktorn
index = 596/532.9;

kostnad_sieve = Cost(d,Param_bottnar(1,:))*kurs*lang*index*bottnar_verklig;
kostnad_valve = Cost(d,Param_bottnar(2,:))*kurs*lang*index*bottnar_verklig;
kostnad_bubble = Cost(d,Param_bottnar(3,:))*kurs*lang*index*bottnar_verklig;

% Skalmassa
Param_skalmassa = [10200 31 0.85
                   12800 73 0.85];

kostnad_colonwall_kol = Cost(mwall_dest(1),Param_skalmassa(1,:))*kurs*lang;
kostnad_colonwall_rostfri = Cost(mwall_dest(2),Param_skalmassa(2,:))*kurs*lang;

Tot_kost_dest = kostnad_sieve+kostnad_colonwall_kol;
% clc
% disp(['Rec    ', num2str(Rec)])
% disp(['Kostnad    ', num2str(Tot_kost_dest)])

%% UTSKRIVNING AV RESULTAT: Separation
disp(' ')
disp(' ')
disp('Destillation:')
disp(' ')
disp('______________________Dimensionering__________________________')
disp(['xf:                               ',num2str(xf, '%.3f')])
disp(['Antal ideala bottnar              ' num2str(bottnar_ideal)])
disp(['Antal verkliga  bottnar           ' num2str(bottnar_verklig)])
disp(['Diameter på tornet (m)            ' num2str(d, '%.2f')])
disp(['Höjd på tornet (m)                ' num2str(h)])
disp(' ')

% Vapour and liquid flowrates in kmol h^-1
disp(['L (förstärkardel) (kmol h-1)      ' num2str(L, '%.2f')])
disp(['V (förstärkardel) (kmol h-1)      ' num2str(V, '%.2f')])
disp(['L (avdrivardel) (kmol h-1)        ' num2str(l, '%.2f')])
disp(['V (avdrivardel) (kmol h-1)        ' num2str(v, '%.2f')])
disp(' ')

% Värmen
disp(['Kondensorvärme (MW)               ' num2str(Q_condensor*10^-6, '%.3f')])
disp(['Återkokarvärme (MW)               ' num2str(Q_reboiler*10^-6, '%.3f')])
disp(' ')

% Factors from correlation
disp(['F_LV                              ' num2str(Flv, '%.3f')])
disp(['C_F (ft s^-1)                     ' num2str(Cf/0.3048)])
disp(' ')

% Kostnad bottnar
disp('______________________Separationkostnader_____________________')
disp(['Kostnad (sieve) (kr)              ' num2str(kostnad_sieve, '%.2f')])
disp(['Kostnad (valve) (kr)              ' num2str(kostnad_valve, '%.2f')])
disp(['Kostnad (bubble) (kr)             ' num2str(kostnad_bubble, '%.2f')])
disp(' ')

% Väggtjocklek och pris på kolonn
disp(['Volym (kolstål) (m^3)             ' num2str(Vwall_dest, '%.2f')])
disp(['Volym (rostfritt stål) (m^3)      ' num2str(Vwall_dest, '%.2f')])
disp(['Massa (kolstål) (kg)              ' num2str(mwall_dest(1), '%.2f')])
disp(['Massa (rostfritt stål) (kg)       ' num2str(mwall_dest(2), '%.2f')])
disp(['Kostnad (kolstål) (kr)            ' num2str(kostnad_colonwall_kol, '%.2f')])
disp(['Kostnad (rostfritt stål) (kr)     ' num2str(kostnad_colonwall_rostfri, '%.2f')])
disp(['Kostnad kolstål + silbottnar (kr) ' num2str(Tot_kost_dest, '%.2f')])


%% Funktioner

% Reaktorfunktioner
function dUdV = PBR_ode(~, U, dHr0,A,E,K1,K2,P,Cp)

R = 8.31447;               % gaskon.

% A = isobutAn (C4H10)
% B = isobutEn (C4H8)
% C = Vätgas   (H2)
% D = vatten   (H2O)

% Hämtar molflöden och temperatur
FA = U(1); 
FB = U(2); 
FC = U(3); 
FD = U(4); 
T = U(5);
F_tot = FA + FB + FC + FD;

% Cp-funktioner beroende av T
CpA = @(T) Cp_calc(T,Cp(1,:));  %Cp för ämne A (isobutan) vid T
CpB = @(T) Cp_calc(T,Cp(2,:));  %Cp för ämne B (isobuten) vid T
CpC = @(T) Cp_calc(T,Cp(3,:));  %Cp för ämne C (H2) vid T
CpD = @(T) Cp_calc(T,Cp(4,:));  %Cp för ämne D (H2O) vid T

%Beräknar partialtryck
PA = P*(FA/F_tot);  %Partialtryck för ämne A (isobutan) för inflöde till reaktor
PB = P*(FB/F_tot);  %Partialtryck för ämne B (isobuten)) för inflöde till reaktor
PC = P*(FC/F_tot);  %Partialtryck för ämne C (H2)) för inflöde till reaktor
% PD = P*(FD/F_tot);  %Partialtryck för ämne D (H2O) för inflöde till reaktor

% Ke som funktion av T
Ke = (2.1*10^7) * exp(-122*10^3/(8.314*T)); % [bar]

% Reaktionskonstant k beroende av T
k = A*exp(-E/(R*T));

% Beräknar reaktionshastigheten r
r = k*(PA - ((PB*PC)/Ke))/(1+K1*PB*((PC)^0.5)+(K2*PC)^0.5); % [mol kg cat.^-1 s^-1 bar^-0.5]

% Massbalanser för A -> B + C
dFA = -r;
dFB =  r;
dFC =  r;
dFD =  0;

% deltaH för A -> B + C
dHr = dHr0 + integral(@(T) CpA(T)-CpB(T)-CpC(T),298.15,T);  % Standardreaktionsentalpin beräknas vid 298.15 K

% Temperaturbalans
dTdV = r*(-dHr) / ((CpA(T)*FA + CpB(T)*FB + CpC(T)*FC + CpD(T)*FD));

dUdV = [dFA 
        dFB 
        dFC 
        dFD 
        dTdV];
end

% Beräkning av Cp vid T med hjälp av coeff.
function Cp = Cp_calc(T, CPcoeff)
    A = CPcoeff(1);
    B = CPcoeff(2);
    C = CPcoeff(3);
    D = CPcoeff(4);

    Cp = A + B*T + C*T.^2 + D*T.^3;

end

% % Beräkning av NTU från epsilon
% function diff = eps(epsilon,Cmm,NTU)
% 
% VL=epsilon;
% HL=(1-exp(-NTU.*(1-Cmm)))./(1-Cmm.*exp(-NTU.*(1-Cmm)));
% 
% diff=HL-VL;
% 
% end

% Separationsfunktioner
function [gamma1, gamma2] = wilson(x1,W12,W21)
    x2 = 1 - x1;

    gamma1 = exp(-log(x1+W12.*x2) + x2.*(W12./(x1+W12.*x2) - W21./(W21.*x1+x2)));
    gamma2= exp(-log(x2+W21.*x1) - x1.*(W12./(x1+W12.*x2) - W21./(W21.*x1+x2)));

end

function P_sat = antoine(T,Ant)

    A = Ant(1);
    B = Ant(2); 
    C = Ant(3);
    P_sat = exp(A-(B./(T+C))); 

end

function res = find_Tb(T,x1,gamma1,gamma2,Ant1,Ant2,P)

A1 = Ant1(1); B1 = Ant1(2); C1 = Ant1(3);
A2 = Ant2(1); B2 = Ant2(2); C2 = Ant2(3);

P1 = exp(A1-(B1./(T+C1)));
P2 = exp(A2-(B2./(T+C2)));

y1 = (x1.*gamma1.*P1)./P;
y2 = ((1-x1).*gamma2.*P2)./P;

res=y1+y2-1;

end

% Kostnadsfunktioner
function kostnad = Cost(S,Param)
    a = Param(1);
    b = Param(2);
    n = Param(3);
    
    kostnad = a + b*S.^n;
end