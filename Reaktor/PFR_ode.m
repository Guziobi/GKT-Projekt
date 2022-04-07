function dUdV = PFR_ode(V, U, Cp, dH0, A, Ea)

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

% Data
R  = 8.31447;               % gaskon.
P  = 1;                     % bar 
k  = A*exp(-Ea./(R.*T));    % mol/kg cat.*s*bar vid 550 C
K1 = 22.9;                  % bar^-1.5
K2 = 7.56;                  % bar^-1

% Hämtar Cp
Cp_A = Cp(1); 
Cp_B = Cp(2); 
Cp_C = Cp(3); 
Cp_D = Cp(4);

%Beräknar partialtryck
PA = P*(FA/F_tot); 
PB = P*(FB/F_tot);
PC = P*(FA/F_tot); 
PD = P*(FD/F_tot);

% Arrenius

% Ke som funktion av T
Ke = (2.1*10^7) * exp(-122000/(R*T)); % bar

% Beräknar rA
rA = (k.* (PA - (PB*PC/Ke))) ./ (1 + (K1*PB*sqrt(PC)) + (sqrt(K2*PC)));

sumFCp = FA*Cp_A + FB*Cp_B + FC*Cp_C + FD*Cp_D;

% Massbalanser för A -> B + C
dFA = -rA;
dFB =  rA;
dFC =  rA;
dFD =  0;

% deltaH för A -> B + C
dH  = dH0 + (Cp_B+Cp_C-Cp_B); 

% deltaH för en viss temp ska beräknas via Cp vid en viss temp här för att
% få en dTdV


% Temperaturbalans
dTdV = rA * -dH / (sumFCp);

dUdV = [dFA 
        dFB 
        dFC 
        dFD 
        dTdV];
end