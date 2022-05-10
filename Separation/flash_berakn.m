p = 101325; % Pa
R = 8.3145;
T = 348;

%Antoinekonstanter A  B  C
Ant1 =  [15.7564 2132.42 -33.15];  % buten
Ant2 =  [15.6782 2154.90 -34.42];  % butan 
Ant3 =  [13.6333 164.90 3.19];  % vätgas
Ant4 =  [18.3036 3816.44 -46.13];  % vatten
%Wilsonfaktorer
W12 = 0.48584; 
W21 = 1.64637; 

% Molmassor
MA = 56.1063;      % g mol-1  Butan
MB = 58.1222;      % g mol-1  Buten
MC = 2.0160;       % g mol-1  Vätgas
MD = 18.0160;      % g mol-1  Vatten

rhoA = 559.0;       % kgm-3 Buten
rhoB = 556.62;      % kgm-3 Butan
rhoD = 974.93;      % kg m-3 Vatten flytande


