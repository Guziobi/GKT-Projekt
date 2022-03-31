%Funktion för beräkning av Cp vid temperaturen T för isobuten (g)

function cp_B=cp_B(T)
A_B=16.05;
B_B=0.2804;
C_B=-1.091e-4;
D_B=9.098e-9;

cp_B=A_B+B_B.*T+C_B.*T.^2+D_B.*T.^3;

end