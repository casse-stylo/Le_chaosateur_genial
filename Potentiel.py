import numpy as np
import matplotlib.pyplot as plt
from math import*



def Kepler(x, y):
    return -1/np.sqrt(x**2 + y**2)
def Henon_Heiles (x,y):
    return 0.5*(x**2 + y**2+2*x**2*y-2/3*y**3)

def Fx (Pot, x, y, h= 1e-3):


    return 1/(12*h) * (-Pot(x+2*h, y)+ 8*Pot(x+h, y) - 8*Pot(x-h, y) + Pot(x-2*h, y))


def Fy (Pot, x, y, h = 1e-3):
    return 1/(12*h) * (-Pot(x,y+2*h)+ 8*Pot(x,y+h) - 8*Pot(x,y-h) + Pot(x,y-2*h))

def f(wn, pot=Kepler):
    return np.array([wn[2], wn[3], -Fx(pot,wn[0], wn[1]), -Fy(pot,wn[0], wn[1])])
