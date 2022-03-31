%Funktion för beräkning av Cp vid temperaturen T för vätgas

function cp_C=cp_C(T)
t=T./1000;

A_C=33.066178;
B_C=-11.363417;
C_C=11.432816;
D_C=-2.772874;
E_C=-1.158558;

cp_C=A_C+B_C.*t+C_C.*t.^2+D_C.*t.^3+E_C.*t^(-2);

end






