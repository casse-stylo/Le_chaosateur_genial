import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*

#test : cercle

def f(wn, pot=Kepler):
    return np.array([wn[2], wn[3], -Fx(pot,wn[0], wn[1]), -Fy(pot,wn[0], wn[1])])

def RK2 (wn, f, h):

    wn12 = wn + h*0.5*f(wn)
    wn12p = f(wn12)

    return wn + h*wn12p

wn = np.array([1,0,0,1])

for i in range(1000) : 
    wn = RK2(wn, f, 1e-2)
    plt.scatter(wn[0],wn[1])

plt.show()