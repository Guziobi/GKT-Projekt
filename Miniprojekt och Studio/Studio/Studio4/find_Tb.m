function res = find_Tb(x1,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P,T)
%Use T,x1,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P
%to calculate y1 and y2



P0_1 = 10.^(A1 - (B1)./(T + C1));
P0_2 = 10.^(A2 - (B2)./(T + C2));

y1 = (P0_1.* x1.* gamma1)./P;
y2 = (P0_2.*(1-x1).*gamma2)./P;

res = y1+y2-1;


% res = [ diff
%         P0_1
%         P0_2
%         y1
%         y2
%         ];

end