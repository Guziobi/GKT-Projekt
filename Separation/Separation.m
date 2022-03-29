% Separation
clc, clear

% Data - separation 
F = 100;           % kmol h-1
P = 760;           % mmHg (1 atm)
xf = 0.50;         % molbråk buten
xd = 0.95;         % destillatbråk 
xb = 0.05;         % bottenbråk
R = ;              % återflödesförhållande
M1 = 56.1063;      % g mol-1
M2 = 58.1222;      % g mol-1

         % A      B       C
Ant1 =  [15.7564 2132.42 -33.15];  % buten
Ant2 =  [15.6782 2154.90 -34.42];  % butan 

W12 = 0.48584; 
W21 = 1.64637; 

rho_butenL = 1;          % kgm-3
rho_butanL = 1;          % kgm-3

liqsurften = 1;         % dyn cm-1

