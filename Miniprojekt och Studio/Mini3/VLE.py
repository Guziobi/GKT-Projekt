from math import *

# Beräknar ångtryck med Antoines lag
def antoine(ABC, T):
    #ABC = [A, B, C]
    P0_i = 10**(ABC[0] - (ABC[1])/(T + ABC[2]))

    return P0_i


# Interpolerar mellan två gissnigar och genererar nästa värde
def interpolation(T_n, T_n1, P, P_n, P_n1):
    T_n2 = T_n1 + ((P - P_n1) * (T_n - T_n1))/(P_n - P_n1)

    return T_n2


# Beräknar aktivitetsfaktorer för de två ämnena
def willson(x1, Lambda12, Lambda21):
    x2 = 1 - x1

    ln_gamma1 = - log(x1 + Lambda12*x2) + x2*(Lambda12 /(x1 + Lambda12*x2) - Lambda21/(Lambda21*x1 + x2))
    ln_gamma2 = - log(x2 + Lambda21*x1) - x1*(Lambda12 /(x1 + Lambda12*x2) - Lambda21/(Lambda21*x1 + x2))

    gamma1 = exp(ln_gamma1)
    gamma2 = exp(ln_gamma2)

    return gamma1, gamma2


def relflyktighet(x1, y1):
    x2 = 1-x1
    y2 = 1-y1
    alpha12 = (y1/x1)/(y2/x2)

    return alpha12


def bubbeltemp(x1, tol, P, Tb1, Tb2, W12, W21, ABC1, ABC2):
   
    diff = 1

    # Beräknar ångtryck vid lättflyktig komponents kokpunkt som första gissning
    T_n = Tb1
    P0_1 = antoine(ABC1, T_n)
    P0_2 = antoine(ABC2, T_n)
    P_n = P0_1*x1 + P0_2*(1-x1)

    # Beräknar ångtryck vid tung komponents kokpunkt som andra gissning
    T_n1 = Tb2
    P0_1 = antoine(ABC1, T_n1)
    P0_2 = antoine(ABC2, T_n1)
    P_n1 = P0_1*x1 + P0_2*(1-x1)

    while diff >= tol:
        # for i in range(10):
        T_n2 = interpolation(T_n, T_n1, P, P_n, P_n1)

        T_n = T_n1
        T_n1 = T_n2

        P0_1 = antoine(ABC1, T_n1)  # Beräkna ångtryck för komponent 1
        P0_2 = antoine(ABC2, T_n1)  # Beräkna ångtryck för komponent 2

        gamma1, gamma2 = willson(x1, W12, W21)

        y1 = (P0_1 * x1 * gamma1)/P  # Beräknar ångfraktionen
        y2 = (P0_2 * (1 - x1) * gamma2)/P

        # Bryter loopen om 1 - sum av y är mindre än tolleransen
        diff = abs((y1+y2) - 1)

        P_n = P_n1
        P_n1 = P0_1*x1*gamma1 + P0_2*(1-x1)*gamma2

    return T_n1, x1, y1



#Flashberäkningar
def rahsfordRice(psi, zi, Ki):
    rice = (zi*(1 - Ki)) / (1 + psi*(Ki - 1))

    return rice

def rahsfordRicePrime(psi, zi, Ki):
    ricePrime = (zi*(1-Ki)**2) / ( (1 + psi*(Ki-1))**2 )

    return ricePrime

def riceRecursion(psi, rice, ricePrime):
    psi2 = psi - (rice/ricePrime)

    return psi2


def flash(z1, P_out, T_out, ABC1, ABC2, tol=0.0001, psi=0.5):
    diff = 1
    z2 = 1-z1

    P01 = antoine(ABC1, T_out)
    P02 = antoine(ABC2, T_out)
    K1 = P01/P_out
    K2 = P02/P_out
    print(K1, K2, '\n')

    for i in range(2):
    #while diff >= tol:

        rice = rahsfordRice(psi, z1, K1) + rahsfordRice(psi, z2, K2)
        print(i, 'rice: ', rice)

        ricePrime = rahsfordRicePrime(psi, z1, K1) + rahsfordRicePrime(psi, z2, K2)
        print('ricePrime: ', ricePrime)
        
        psi2 = riceRecursion(psi, rice, ricePrime)
        
        diff = abs((psi2-psi)/psi)
        psi = psi2
        print('psi: ', psi)
        print('diff: ', diff, '\n')

    print('klar')
    

    return psi, diff


def flash2(F, L, z1, P_out, T_out, ABC1, ABC2):
    z2 = 1-z1

    P01 = antoine(ABC1, T_out)
    P02 = antoine(ABC2, T_out)
    x1 = (F*z1)/ ((F-L)*(P01/P_out) + L)
    x2 = (F*z2)/ ((F-L)*(P02/P_out) + L)
    xsum = x1 + x2

    return x1, x2, xsum