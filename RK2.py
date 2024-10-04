import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*

#test : cercle

def f(wn, pot=Kepler):
    return np.array([wn[2], wn[3], -Fx(pot,wn[0], wn[1]), -Fy(pot,wn[0], wn[1])])

def RK2 (wn, f, h, pot):

    wn12 = wn + h*0.5*f(wn,pot)
    wn12p = f(wn12,pot)

    return wn + h*wn12p

