%Funktion för beräkning av Cp vid temperaturen T för flytande vatten

function cp_Dl=cp_Dl(T)
t=T./1000;

A_Dl=-203.606;
B_Dl=1523.29;
C_Dl=-3196.413;
D_Dl=2474.445;
E_Dl=3.855326;

cp_Dl=A_Dl+B_Dl.*t+C_Dl.*t.^2+D_Dl.*t.^3+E_Dl.*t^(-2);

end