clear, clc
Tf = 70; % C
P=760.0021; %mmHg
xfa = 0.45; % konc. av A i vätskeinflödet
zfa = xfa;   % toaltkonc av A i inflödet
xd=0.8;
xw=0.1;     %komponent x i nedre delen av kolonnen (W och B motsvarar samma del)
F=100;                  %feed
R=2;                    %återflödesförhållande


%Antoine konstanter
Ant_a=[18.9119 3803.93 231.47];
Ant_b=[17.5439 3166.34 193.00];

%Enatalpi och cp
Hvap_flow = 41000; %J mol-1
Cp_flow = 130; %J mol-1 K-1

%q-värde
Tguess = 10;

TB_flow = fsolve(@(T)find_Tb(xfa,1,1,Ant_a(1),Ant_a(2),Ant_a(3),Ant_b(1),Ant_b(2),Ant_b(3),P,T), Tguess);

q = (Hvap_flow + Cp_flow*(TB_flow - Tf))/Hvap_flow;

% totalbalans + Komponentbalans för att beräkna strömarnas storlekar
D = (F.*zfa - F.*xw)./(xd - xw);
W = F-D; %Samma som B i powerpoint-slidesen

% balances over condenser.
L = R.*D;
V = D.*(R + 1);

l = L + q.*F; % L-streck
v = V - (1 - q).*F; % V-streck

% relativ volatilitet
alpha = 2.08;

%återkokare
y0= (xw*alpha)/(1 + xw*(alpha-1));
x11= (v.*y0 + W.*xw)./l;         % KB återkokare

% Avdrivardel Stegning till feedbotten:
x(1)=x11;
i=0;


while x<xfa
    i=i+1;
    y(i)=(x(i)*alpha)/(1 + x(i)*(alpha-1));
    x(i+1)= (v.*y(i) + W.*xw)./l; % Komponentbalans för avdrivardelen
end

m=i+1;

%Förstärkare
while y(i)<xd
    x(i+1) = (V.*y(i) + W.*xw - F.*zfa)./L;   %Komponentbalans över förstärkardelen
    i=i+1;
    y(i)=(x(i)*alpha)/(1 + x(i)*(alpha-1)); %Jv samband
end

%% Bottendiameter

% Sammansättning ut ur återkokare
v_x1 = y0;              % Komponent A
v_x2 = 1-y0;            % Komponent B

% Sammansättning in i återkokare
l_x1 = x(1);
l_x2 = 1 - x(1);

% Konstanter
%Densiteter
rho_e = 772; %kg m-3
rho_p = 796; %kg m-3
M_e = 46.06844e-3; %kg/mol
M_p = 60.0952e-3; %kg/mol

% Ändra detta
M_L = l_x1*M_e + l_x2*M_p;
M_V = v_x1*M_e + v_x2*M_p;

rho_L = ((l_x1*M_e)/(l_x1*M_e + l_x2*M_p))*rho_e + ((l_x2*M_p)/(l_x1*M_e + l_x2*M_p))*rho_p; % kg m-3
rho_V = 1.8; %kg m-3

% Flödesfaktorer
surftention = 24; %dyn cm-1
Fst = (surftention/20)^0.2;
vaporveloc = 0.7;
Ff = 1; %nono-foaming
Fha = 1; %Hålen är bra


Flv = ((l*M_L)/(v*M_V)) * sqrt(rho_V/rho_L);

% från diagramet
Cf = 0.38*0.3048; %ft/s -> m/s

%Flödningsparametern
C = Fst*Ff*Fha*Cf;

% Ånghastigheten vid flödning
Uf = C*sqrt((rho_L-rho_V)/rho_V);

%Ånghastighet
U = vaporveloc*Uf;

%Aktiv area
Aaktiv = (V*(1/3.6)*M_V*(1/rho_V))/U;

%Total area
Atot = Aaktiv/0.8;

%Bottendiameter
d = 2*sqrt(Atot/pi);

%% Andra botten
TB_2= fsolve(@(T)find_Tb(x(2),1,1,Ant_a(1),Ant_a(2),Ant_a(3),Ant_b(1),Ant_b(2),Ant_b(3),P,T), Tguess);

%% Värmebehov
% återkokare
qw= v.*Hvap_flow*(1/60)*(1/60);        % Ångbildningsvärme i kW
% kondensor
qc= V.*Hvap_flow*(1/60)*(1/60);        % Kylning i kW


%% Display
disp(['bottnar nedre del    ' num2str(m)])
disp(['bottnar övre del     ' num2str(i-m)])

disp('')
disp(['bottnar totalt       ' num2str(i)])
disp(['Bottendiameter       ' num2str(d)])

disp(' ')
disp(['Förstärkardel L     ' num2str(L)])
disp(['Förstärkardel V     ' num2str(V)])
disp(['Avdrivardel L       ' num2str(l)])
disp(['Avdrivardel V       ' num2str(v)])

disp(' ')
disp(['q-värde     ' num2str(q)])
disp(['Värme kondensator     ' num2str(qc*1e-3), ' MW'])
disp(['Värme återkokare      ' num2str(qw*1e-3), ' MW'])

disp(' ')
disp(['Flv (x10)     ' num2str(Flv*10)])
disp(['Cf     ' num2str(Cf/0.3048)])

disp(' ')
disp(['Vätska, x1 vid andra botten    ' num2str(x(2))])
disp(['Ånga, y1 vid andra botten    ' num2str(y(2))])
disp(['Temperatur vid andra botten    ' num2str(TB_2)])





figure(1);
plot(1:i,y,'r');
hold on
plot(1:i,x);
legend('y1','x1');
ylabel('x1,y1');
xlabel('Botten nr (Räknat från återkokaren)');

% Minsta återflödesförhållandet

yfa=5/3-5/(3+4.5*xfa);
Rmin=(xd-yfa)/(yfa-xfa);            % (LV)min = Rmin/(Rmin+1)
disp(['Minsta återflödesförhålland:    ',num2str(Rmin)])