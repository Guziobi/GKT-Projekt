clear, clc;

% Heat exchanger data
L = 3.75;       % [m] tube lenght
n = 68;         % [#] number of tubes
Di = 0.088;     % [m] inner tube diameter
Do = 0.094;     % [m] outer tube diameter
kwall = 35;     % [W m-1 K-1] tube lenght

% Stream data
mprc = 3.75;    % [kg s-1] Mass flow
mprh = 2.0625;  % [kg s-1] Mass flow
TCin = 15;      % [C] Temp in
THin = 90;      % [C] Temp in
hc = 300;       % [W m-2 K-1] Heat transfer coefficients
hh = 700;       % [W m-2 K-1] Heat transfer coefficients
cpc = 1900;     % [J kg-1 K-1] Specific heat
cph = 4200;     % [J kg-1 K-1] Specific heat


%% Case 1

ri = Di/2;
ro = Do/2;
Ac = pi*Di*L*n;
Ah = pi*Do*L*n;

C = [mprc*cpc mprh*cph];
Cmin = min(C);
Cmax = max(C);

UA = ( (hh*Ah)^-1 + (hc*Ac)^-1 + log(ro/ri)/(2*pi*kwall*L*n) )^-1;
% epsilon = @(NTU)(1 - exp(-NTU.*(1 - Cmin/Cmax)))/(1 - (Cmin/Cmax).*exp(-NTU.*(1 - Cmin/Cmax)));

NTU = UA/Cmin;
Ckvot = Cmin/Cmax;

e = 0.61;

q = e*Cmin*(THin - TCin);

TCut = q/(mprc*cpc) + TCin;
THut = THin - q/(mprh*cph);

% Display results
disp('Case 1')
disp(['UA           ' num2str(UA) ' W/K'])
disp(['Cmin         ' num2str(Cmin) ' W/K'])
disp(['Ckvot        ' num2str(Ckvot)])
disp(['NTU          ' num2str(NTU)])
disp(['Epsilon      ' num2str(e)])
disp(['TCout        ' num2str(TCut) ' C'])
disp(['THout        ' num2str(THut) ' C'])
disp(' ')

%% Case 2

UA = ( (hh*Ah*1.5)^-1 + (hc*Ac*1.5)^-1 + log(ro/ri)/(2*pi*kwall*L*n) )^-1;
NTU = UA/Cmin;
e = 0.65;

q = e*Cmin*(THin - TCin);
TCut = q/(mprc*cpc) + TCin;
THut = THin - q/(mprh*cph);

% Display results
disp('Case 2')
disp(['UA           ' num2str(UA) ' W/K'])
disp(['Cmin         ' num2str(Cmin) ' W/K'])
disp(['Ckvot        ' num2str(Ckvot)])
disp(['NTU          ' num2str(NTU)])
disp(['Epsilon      ' num2str(e)])
disp(['TCout        ' num2str(TCut) ' C'])
disp(['THout        ' num2str(THut) ' C'])
disp(' ')

%% Case 3

C = [mprc*1.6*cpc mprh*cph];
Cmin = min(C);
Cmax = max(C);

Ckvot = Cmin/Cmax;

UA = ( (hh*Ah)^-1 + (hc*Ac)^-1 + log(ro/ri)/(2*pi*kwall*L*n) )^-1;
NTU = UA/Cmin;
e = 0.6;

q = e*Cmin*(THin - TCin);
TCut = q/(mprc*1.6*cpc) + TCin;
THut = THin - q/(mprh*cph);

% Display results
disp('Case 3')
disp(['UA           ' num2str(UA) ' W/K'])
disp(['Cmin         ' num2str(Cmin) ' W/K'])
disp(['Ckvot        ' num2str(Ckvot)])
disp(['NTU          ' num2str(NTU)])
disp(['Epsilon      ' num2str(e)])
disp(['TCout        ' num2str(TCut) ' C'])
disp(['THout        ' num2str(THut) ' C'])


