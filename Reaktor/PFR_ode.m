function dUdV = PFR_ode(W, U, dHr0,k,K1,K2,P,Cp)

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
PD = P*(FD/F_tot);  %Partialtryck för ämne D (H2O) för inflöde till reaktor

% Ke som funktion av T
Ke = (2.1*10^7) * exp(-122*10^3/(8.314*T)); % [bar]

% Beräknar r
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