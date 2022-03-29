% Program to solve a system of coupled ordinary differential equations
% 
% Reaction carried out in Batch Reactor
% Reaction:        ONCB + 2NH3 -> NA + ACl
% Rate Epression:  -r_ONCB = k C_ONCB C_NH3
% 
clear; clc;

%Set values of constants
R = 8.31446261815324;   % m3 Pa K−1 mol−1
Ea = 47.17e3;           % J mol-1
DHr = -2.47e6;          % J mol-1
A = 37.61;              % dm3 mol-1 min-1
UA = 150.0e3;           % J min-1 K-1
V = 5.12e3;             % [dm3] volume of reactor
Cp_ONCB = 167.4;        % J mol-1 K-1
Cp_H2O = 75.3;          % J mol-1 K-1
Cp_NH3 = 35.1;          % J mol-1 K-1
T_a = 25 + 273;         % [K]


%Init numbers of mole
N0_ONCB = 9.04e3;  N0_NH3 = 33.0e3; N0_H2O = 103.7e3;   % mol
N0_NA = 0; N0_ACl = 0;                                  % mol
T0 = 175+273;                                           % [K] init temp
U0 = [N0_ONCB N0_NH3 N0_NA N0_ACl T0];                  % Init values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start and Stop times
t_start=0;              % start value for t
t_final=24*60;          % final value for t [min]
tspan = [t_start t_final];

%[t_out,U_out] = ode45(@ode_function,[t_range],[initial_values],[options],constants)
[t, U]=ode15s(@ode_accident_cooling, [t_start t_final], U0,[], R, Ea, DHr, A, UA, V, N0_H2O, Cp_ONCB, Cp_H2O, Cp_NH3);

%extract N_ONCB etc....
N_ONCB=U(:,1); N_NH3=U(:,2); N_NA=U(:,3); N_ACl=U(:,4); T=U(:,5);

%Concentrations
C_ONCB = N_ONCB./V;
C_NH3  = N_NH3./V;

%rate
k = A*exp(-Ea/R./T);
r_ONCB = -k.*C_ONCB.*C_NH3;

%calculate conversion of A
XA=(N0_ONCB - N_ONCB)./N0_ONCB;

%calculate Q's
Qmax_r = -UA.*(T_a - T);
Q_g = r_ONCB.*V.*DHr;

%plot conversion of A versus time
figure(1)
plot(t,XA)
xlabel('t [min]')
ylabel('Conversion of A')

%plot T versus time
figure(2)
plot(t,(T-273))
xlabel('t [min]')
ylabel('T [C^{o}]')

%plot T versus time
figure(3)
plot(t,Qmax_r*10^-3, t,Q_g*10^-3)
xlabel('t [min]')
ylabel('Q [kJ min^-1]')

