function dUdV = ode_pfr(U, R, Ea, Cp, T)

NA = U(1); NB = U(2); NC = U(3); ND = U(4);
Cp_A = Cp(1); Cp_B = Cp(2); Cp_C = Cp(3); Cp_D = Cp(4);

sumNCp = NA*Cp_A + NB*Cp_B + NC*Cp_C + ND*Cp_D;


PA = NA * R * T / V;
PB = NB * R * T / V;
PC = NC * R * T / V;

rA = ( k* (PA - (PB*PC/Ke))) / (1 + K1*PB*sqrt(PC) + sqrt(K2*PC));

% Massbalanser

dNA = -rA * dV;
dNB =  rA * dV;
dNC =  rA * dV;
dND = 0;


% Temperaturbalans
dTdV = r * -dH / (sumNCp)

dUdV = [dNA dNB dNC dND dTdV];
end