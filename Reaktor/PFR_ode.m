function dUdV = PFR_ode(U, Cp)

% A = isobutAn (C4H10)
% B = isobutEn (C4H8)
% C = Vätgas   (H2)
% D = vatten   (H2O)

% Data
R  = 8.31447; % gaskon.
P  = 1;       % bar 
k  = 0.0596;  % mol/kg cat.*s*bar vid 550 C
K1 = 22.9;    % bar^-1.5
K2 = 7.56;    % bar^-1

% Hämtar molflöden och temperatur
FA = U(1); FB = U(2); 
FC = U(3); FD = U(4); 
T = U(5);
F_tot = FA + FB + FC + FD;

% Hämtar Cp
Cp_A = Cp(1); Cp_B = Cp(2); 
Cp_C = Cp(3); Cp_D = Cp(4);

%Beräknar partialtryck
PA = P*(FA/F_tot); PB = P*(FB/F_tot);
PC = P*(FA/F_tot); PD = P*(FD/F_tot);

% Ke som funktion av T
Ke = @(T)(2.1*10^7) * exp(-122000/(R*T)) % bar


% Beräknar rA
rA = (k* (PA - (PB*PC/Ke(T)))) / (1 + K1*PB*sqrt(PC) + sqrt(K2*PC));

sumFCp = FA*Cp_A + FB*Cp_B + FC*Cp_C + FD*Cp_D;

% Massbalanser för A -> B + C
dFA = -rA;
dFB =  rA;
dFC =  rA;
dFD =  0;

% deltaH för en viss temp ska beräknas via Cp vid en viss temp här för att
% få en dTdV


% Temperaturbalans
dTdV = r * -dH / (sumFCp)

dUdV = [dFA dFB dFC dFD dTdV];
end