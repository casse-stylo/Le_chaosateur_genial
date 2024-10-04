import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*


def RK4 (wn, f, h, pot):
    
    k1 = f(wn, pot)
    k2 = f(wn + h/2 * k1, pot)
    k3 = f(wn+ h/2 * k2, pot)
    k4 = f(wn+h*k3, pot)



    return wn + h/6 * (k1 + 2*k2 + 2*k3 + k4)

