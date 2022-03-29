function y1 = flash_ode(T, P, z1, A1, B1, C1, W12, W21)

% F = V + L
% (F*z1 - L*x1)./V = y1 
% V = 0.5.*F


P0_1 = 10.^(A1 - (B1)./(T + C1));
% P0_2 = 10.^(A2 - (B2)./(T + C2));

x1 = (z1)./(0.5 + (0.5.*P0_1)/P);
x2 = 1 - x1;

gamma1 = exp(- log(x1 + W12.*x2) + x2.*(W12./(x1 + W12.*x2) - W21./(W21.*x1 + x2)));

y1_start = (z1 - 0.5.*x1)./0.5;
y1_jmk = (P0_1.* x1.* gamma1)./P;

y1 = [  y1_start
        y1_jmk];


end
