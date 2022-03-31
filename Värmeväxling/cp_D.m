%Funktion för beräkning av Cp vid temperaturen T för vattenånga

function cp_D=cp_D(T)
t=T./1000;

A_D=30.092;
B_D=6.832514;
C_D=6.793435;
D_D=-2.53448;
E_D=0.082139;

cp_D=A_D+B_D.*t+C_D.*t.^2+D_D.*t.^3+E_D.*t^(-2);

end