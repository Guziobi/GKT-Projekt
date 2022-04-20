% KAA146 Grundläggnde Kemiteknik, Projektarbete Grupp 4
% Modell för reaktorerna

%% Massbalans över reaktorn
clc
clear

% Fx_x är molfflöden i kmol/h
% Mx_x är massflöde i ton/h

M_A = 58.124;           %Molmassa [kmol/kg] Butan C4H10
M_B = 56.1063;          %Molmassa [kmol/kg] Buten C4H8
F11_onskad = 6250./M_B; %Önskad produktion [kmol/h]
F1_A = 10000./M_A;      %Basis [kmol/h]
F3_B = (F1_A.*0.45);     
F3_A = (F1_A-F3_B);
F5_B = (F3_A.*0.45 + F3_B);
F5_A = (F1_A-F5_B);
% SEP molflöden
F10_A = (0.95.*F5_A);
F10_B = (0.05.*F5_B);
F12_B = (0.05.*F10_B);
F12_A = (0.05.*F10_A);
F13_B = (0.95.*F10_B);
F13_A = (0.95.*F10_A);
F0_A = (F1_A-F13_A);
F1_B = (F13_B);
F11_A = (0.05.*F5_A);
F11_B = (0.95.*F5_B);


kvot = (F11_onskad./F11_B);

% Omräkning från bas till faktiskt flöde med kvoten
F1 = [F1_A F1_B].*kvot;
F3 = [F3_A F3_B].*kvot;
F5 = [F5_A F5_B].*kvot;
F10 = [F10_A F10_B].*kvot;
F11 = [F11_A F11_B].*kvot;
F12 = [F12_A F12_B].*kvot;
F13 = [F13_A F13_B].*kvot;
F0 = F0_A.*kvot;
