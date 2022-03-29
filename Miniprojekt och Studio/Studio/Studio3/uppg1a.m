%% Konstanter
clc, clear
% 1 = acetone  2 = methanol
%Wilson parameters
W12=0.65675;
W21=0.77204;
%Antoine constants for degC, mmHg, log10
A1=7.02447;  B1=1161.00;  C1=224;  %acetone
A2=7.87863;  B2=1473.11;  C2=230;  %methanol
%total pressure
P=760;  %mmHg

steg = 100-1;
h = 1/steg;

%% Beräkningar

x1=0:h:1;
x2=1-x1;

Tstart= linspace(56.05,78.37,(steg+1));  %temperature at which to start the search

%activity coefficients at x1
gamma1=exp(- log(x1 + W12.*x2) + x2.*(W12./(x1 + W12.*x2) - W21./(W21.*x1 + x2)));
gamma2=exp(- log(x2 + W21.*x1) - x1.*(W12./(x1 + W12.*x2) - W21./(W21.*x1 + x2)));


%use fsolve function to find bubble point temperature (Tb) for x1
%find_Tb is a function we need to create that will check if a certain value of T satisfies y1+y2-1=0 

%current value of x1 and other constants are passed to find_Tb
Tb=fsolve(@(T)find_Tb(T,x1,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart);


P1 = 10.^(A1 - (B1)./(Tb + C1));
y1 = (P1.*gamma1.*x1)./P;

figure(1);
plot(x1,Tb, y1,Tb)
axis([0 1 55 65])

figure(2);
plot(x1,y1, x1,x1, 'k')
axis([0 1 0 1])