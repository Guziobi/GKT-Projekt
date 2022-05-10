clc, clear
P = 101325; % Pa
R = 8.3145;
T = 348;

%Antoinekonstanter A  B  C
Ant1 =  [15.7564 2132.42 -33.15];  % buten
Ant2 =  [15.6782 2154.90 -34.42];  % butan 
Ant3 =  [13.6333 164.90 3.19];  % vätgas
Ant4 =  [18.3036 3816.44 -46.13];  % vatten


tol = 0.01;
xff1 = 0.35; % Molbråk i feed till flash 1
xff1 = 0.35; % Molbråk i feed till flash 2

%Ungefärliga kokpunkter används som gissningar
Tb_1 = 273.15+20; % Buten @ 6atm [K]
Tb_2 = 272.7; % Butan @ 1atm [K]
Tb_3 = 28.3; % Vätgas @ 6atm [K]
Tb_4 = 273.15+100; % Vatten @ 1atm [K]

% Molmassor
MA = 56.1063;      % g mol-1  Butan
MB = 58.1222;      % g mol-1  Buten
MC = 2.0160;       % g mol-1  Vätgas
MD = 18.0160;      % g mol-1  Vatten

rhoA = 559.0;       % kgm-3 Buten
rhoB = 556.62;      % kgm-3 Butan
rhoD = 974.93;      % kg m-3 Vatten flytande


X_flash1 = linspace(0.1,0.9,100);
[T_bp, x, y] = bubbeltemp(X_flash1, tol, P, Tb_2, Tb_4, Ant2, Ant4);

% % Jämviktkurva och jämviktsplot
% xeq = 0:0.001:1;    
% [gamma1, gamma2] = wilson(xeq,W12,W21);
% Tstart = linspace(-6.3+273.15,-1+273.15,1001);
% TBeq=fsolve(@(T)find_Tb(T,xeq,gamma1,gamma2,Ant1,Ant2,P),Tstart,options);
% Psat1 = antoine(TBeq, Ant1);                                            % Ångtryck
% yeq = (gamma1.*xeq.*Psat1)./P;
% 
% figure(4);
% plot(xeq,xeq)   % Referenslinje
% hold on
% plot(xeq,yeq)
% xlabel('x_1'), ylabel('y_1')
% axis([0 1 0 1])
% 
% legend('Referenslinje', 'Jämviktskurva','Location','northwest')



function P_sat = antoine(T,Ant)

    A = Ant(1);
    B = Ant(2); 
    C = Ant(3);
    P_sat = exp(A-(B./(T+C))); 

end

% Interpolerar mellan två gissnigar och genererar nästa värde
function T_n2 = interpolation(T_n, T_n1, P, P_n, P_n1)
    T_n2 = T_n1 + ((P - P_n1) * (T_n - T_n1))/(P_n - P_n1);
end


function [T_n1, x1, y1] = bubbeltemp(x1, tol, P, Tb1, Tb2, ABC1, ABC2);
   
    diff = 1;

    % Beräknar ångtryck vid lättflyktig komponents kokpunkt som första gissning
    T_n = Tb1;
    P0_1 = antoine(T_n, ABC1);
    P0_2 = antoine(T_n, ABC2);
    P_n = P0_1*x1 + P0_2*(1-x1);

    % Beräknar ångtryck vid tung komponents kokpunkt som andra gissning
    T_n1 = Tb2;
    P0_1 = antoine(T_n, ABC1);
    P0_2 = antoine(T_n, ABC2);
    P_n1 = P0_1*x1 + P0_2*(1-x1);

    for i = 1:10
    %while diff >= tol
        % for i in range(10):
        T_n2 = interpolation(T_n, T_n1, P, P_n, P_n1);

        T_n = T_n1;
        T_n1 = T_n2;

        P0_1 = antoine(T_n, ABC1);  % Beräkna ångtryck för komponent 1
        P0_2 = antoine(T_n, ABC2);  % Beräkna ångtryck för komponent 2

        %gamma1, gamma2 = wilson(x1, W12, W21);
        gamma1 = 1;
        gamma2 = 1;

        y1 = (P0_1 * x1 * gamma1)/P;  % Beräknar ångfraktionen
        y2 = (P0_2 * (1 - x1) * gamma2)/P;

        % Bryter loopen om 1 - sum av y är mindre än tolleransen
        diff = abs((y1+y2) - 1);

        P_n = P_n1;
        P_n1 = P0_1*x1*gamma1 + P0_2*(1-x1)*gamma2;
    end
 end