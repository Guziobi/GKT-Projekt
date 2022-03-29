%% Konstanter
clc, clear
% 1 = acetone  2 = methanol
%Wilson parameters
W12=0.65675;
W21=0.77204;
%Antoine constants for degC, mmHg, log10
A1=7.02447;  B1=1161.00;  C1=224;  %acetone
A2=7.87863;  B2=1473.11;  C2=230;  %methanol
%total pressure
P=760;  %mmHg

Tstart = 50;
z1 = 0.5;


T = fsolve(@(T)flash_ode(T, P, z1, A1, B1, C1, W12, W21),Tstart);

