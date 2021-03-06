clear, clc

P=760; %mmHg
q= 0.8;
zfa=0.35;   % toaltkonc av A i inflödet
xfa=0.3063; % konc. av A i vätskeinflödet
xd=0.9;
xw=0.1;     %komponent x i nedre delen av kolonnen (W och B motsvarar samma del)
F=250;                  %feed
R=3;                    %återflödesförhållande

%Antoine konstanter
Ant_a=[6.90565 1211.033 220.790];
Ant_b=[6.95464 1344.800 219.482];


% totalbalans + Komponentbalans för att beräkna strömarnas storlekar
D = (F.*zfa - F.*xw)./(xd - xw);
W = F-D; %Samma som B i powerpoint-slidesen

% balances over condenser.
L = R.*D;
V = D.*(R + 1);

l = L + q.*F; % L-streck
v = V - (1 - q).*F; % V-streck

%återkokare
y0=5/3-5/(3+4.5*xw);
x11= (v.*y0 + W.*xw)./l;         % KB återkokare

% Avdrivardel Stegning till feedbotten:
x(1)=x11;
i=0;


while x<xfa
    i=i+1;
    y(i)=5/3-5/(3+4.5*x(i));
    x(i+1)= (v.*y(i) + W.*xw)./l; % Komponentbalans för avdrivardelen
end

m=i+1;

%Förstärkare
while y(i)<xd
    x(i+1)= (V.*y(i) + W.*xw - F.*zfa)./L;   %Komponentbalans över förstärkardelen
    i=i+1;
    y(i)=5/3-5/(3+4.5*x(i)); %Jv samband
end


%% värmebehov återkokare
%Beräkning av temperaturen i återkokaren  (Bubbelpunkt)
Tguess = 100;
TK= fsolve(@(T)find_Tb(xw,1,1,Ant_a(1),Ant_a(2),Ant_a(3),Ant_b(1),Ant_b(2),Ant_b(3),P,T), Tguess);


% Sammansättning ut ur återkokare
x1=y0;              % Komponent A
x2=1-y0;            % Komponent B

%Beräkning av entalpier i gasfas
aAg=0.69381e5;
bAg=0.6752e1;
cAg=0.13199;
aBg=0.31596e5;
bBg=0.15841e2;
cBg=0.15429;

HA=aAg+bAg*TK+cAg*TK^2;
HB=aBg+bBg*TK+cBg*TK^2;
Hblandningg = HA.*x1 + HB.*x2;

%Entalpier vätskefas
aAv=0.19534e5;
bAv=0.63711e2;
cAv=0.12206;
aBv=-0.12588e5;
bBv=0.14150e2;
cBv=0.23130;

hA=aAv+bAv*TK+cAv*TK^2;
hB=aBv+bBv*TK+cBv*TK^2;
hblandningv = hA.*x1 + hB.*x2;

qw= v.*(Hblandningg-hblandningv)*(1/60)*(1/60);        % Ångbildningsvärme i kW

disp(['överförd effekt   ' num2str(round(qw)),'  W'])
disp(['bottnar nedre del    ' num2str(m)])
disp(['bottnar övre del     ' num2str(i-m)])
disp(['bottnar totalt       ' num2str(i)])

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