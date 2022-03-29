import numpy as np
import matplotlib.pyplot as plt
from VLE import *

P = 1520  # mmHg
T = 87 # C
Xf = 0.35
tol = 0.001
Tb_m = 64.7
Tb_e = 77.1

A1 = 7.87863
B1 = 1473.11
C1 = 230.0

A2 = 7.09803
B2 = 1238.71
C2 = 217.0

ABC1 = [A1, B1, C1]
ABC2 = [A2, B2, C2]

Lambda12 = 0.62551
Lambda21 = 0.49384


TBP, x1, y1 = bubbeltemp(Xf, tol, P, Tb_m, Tb_e, Lambda12, Lambda21, ABC1, ABC2)
alpha = relflyktighet(x1,y1)
print(' Tryck: {}\n Bubbelpunktstemp: {:.4f} \n Andel tung komponent: {:.4f}\n Andel flyktig komponent: {:.4f}\n Relativ flyktighet: {:.4f}'.format(P, TBP, x1, y1, alpha))

TBP, x1, y1 = bubbeltemp(Xf, tol, (P*2), Tb_m, Tb_e, Lambda12, Lambda21, ABC1, ABC2)
alpha = relflyktighet(x1,y1)
print(' Tryck: {}\n Bubbelpunktstemp: {:.4f} \n Andel tung komponent: {:.4f}\n Andel flyktig komponent: {:.4f}\n Relativ flyktighet: {:.4f}'.format(P*2, TBP, x1, y1, alpha))

x1 = []
y1 = []
TBP = []

X = list(np.linspace(0,1,100))
for xi in X:
    T_bp, x, y = bubbeltemp(xi, tol, P, Tb_m, Tb_e, Lambda12, Lambda21, ABC1, ABC2)
    TBP.append(T_bp)
    x1.append(x)
    y1.append(y)


def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

azeotrop = intersection(x1, y1)
print(azeotrop)

figure, axis = plt.subplots(1, 2)

axis[0].plot(x1, TBP)
axis[0].plot(y1, TBP)
axis[0].legend(['Bubbelpunkt', 'Daggpunkt'])


axis[1].plot(x1, y1, 'red')
axis[1].plot(x1, x1, 'black')
axis[1].legend(['Jämviktskurva'])
plt.show()