function dUdt=ode_accident_cooling(t, U, R, Ea, DHr, A, UA, V, N_H2O, Cp_ONCB, Cp_H2O, Cp_NH3)
% function file containing differential equations
%extract 
N_ONCB=U(1); N_NH3=U(2); N_NA=U(3); N_ACl=U(4); T=U(5);

sumNCp = N_ONCB*Cp_ONCB + N_H2O*Cp_H2O + N_NH3*Cp_NH3;

%calculate concentrations
C_ONCB=N_ONCB/V;  C_NH3=N_NH3/V;

%calculate reaction rate
k = A*exp(-Ea/R/T);
r_ONCB = -k*C_ONCB*C_NH3;

%differential equations
dN_ONCBdt = r_ONCB * V;
dN_NH3dt  = 2 * r_ONCB * V;
dN_NAdt   = -r_ONCB * V;
dN_ACldt  = -r_ONCB * V;
dTdt      = 0;

%Assemble Differential Equations
dUdt=[dN_ONCBdt
      dN_NH3dt
      dN_NAdt
      dN_ACldt
      dTdt
      ];
end