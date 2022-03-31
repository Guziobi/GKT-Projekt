%Funktion för beräkning av Cp vid temperaturen T för isobutan (g)

function cp_A = cp_A(T)
A_A=-1.39;
B_A=0.3847;
C_A=-1.846e-4;
D_A=2.895e-08;

cp_A=A_A+B_A.*T+C_A.*T.^2+D_A.*T.^3;

end