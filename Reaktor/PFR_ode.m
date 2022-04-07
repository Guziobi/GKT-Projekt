function dUdV = PFR_ode(W, U, dH0,k,K1,K2,P)

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

<<<<<<< HEAD
CPcoeff_H2O = [72.43 1.039*10^-2 -1.497*10^-6 0];
CPcoeff_H2 = [27.14 0.009274 -1.38*10^-5 7.645*10^-9];
CPcoeff_ISOBUTAN = [-1.39 0.3847 -1.846*10^-4 2.895*10^-8];
CPcoeff_ISOBUTEN = [16.05 0.2804 -1.091*10^-4 9.098*10^-9];
=======
% Data
R  = 8.31447;               % gaskon.
P  = 1;                     % bar 
k  = A*exp(-Ea./(R.*T));    % mol/kg cat.*s*bar vid 550 C
K1 = 22.9;                  % bar^-1.5
K2 = 7.56;                  % bar^-1
>>>>>>> 96bd975e9251b4d2352214d5dba87c9bb60fd6f5

CPcoeff = CPcoeff_ISOBUTAN;
CpA = @(T) Cp_calc(T,CPcoeff);

CPcoeff = CPcoeff_ISOBUTEN;
CpB = @(T) Cp_calc(T,CPcoeff);

CPcoeff = CPcoeff_H2;
CpC = @(T) Cp_calc(T,CPcoeff);

CPcoeff = CPcoeff_H2O;
CpD = @(T) Cp_calc(T,CPcoeff);

%Beräknar partialtryck
PA = P*(FA/F_tot); 
PB = P*(FB/F_tot);
PC = P*(FC/F_tot); 
PD = P*(FD/F_tot);

% Ke som funktion av T
<<<<<<< HEAD
Ke = (2.1*10^7) * exp(-122*10^3/(8.314*T)); % bar
=======
Ke = (2.1*10^7) * exp(-122000/(R*T)); % bar
>>>>>>> 96bd975e9251b4d2352214d5dba87c9bb60fd6f5

% Beräknar rA
r = k*(PA - ((PB*PC)/Ke))/(1+K1*PB*((PC)^0.5)+(K2*PC)^0.5);

% sumFCp = FA*Cp_A + FB*Cp_B + FC*Cp_C + FD*Cp_D;

% Massbalanser för A -> B + C
dFA = -r;
dFB =  r;
dFC =  r;
dFD =  0;

% deltaH för A -> B + C
dHr = dH0 + integral(@(T) CpA(T)-CpB(T)-CpC(T),298.15,T);

% Temperaturbalans
dTdV = r*(-dHr) / ((CpA(T)*FA + CpB(T)*FB + CpC(T)*FC + CpD(T)*FD));

dUdV = [dFA 
        dFB 
        dFC 
        dFD 
        dTdV];
end