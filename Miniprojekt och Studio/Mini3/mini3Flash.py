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

Xsum = []
#X1 = []
#X2 = []

X = list(np.linspace(0,99,100))
#print(X)
for L in X:
    xsum, x1, x2 = flash2(100, L, Xf, P, T, ABC1, ABC2)
    Xsum.append(xsum)
    #X1.append(x1)
    #X2.append(x2)


# print(Xsum)
plt.plot(X, Xsum)
plt.show()
