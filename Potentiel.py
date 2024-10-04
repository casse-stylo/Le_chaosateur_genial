import numpy as np
import matplotlib.pyplot as plt
from math import*

def Kepler(x, y):
    return -1/sqrt(x**2 + y**2)
def Henon_Heiles (x,y):
    return 0.5*(x**2 + y**2+2*x**2*y-2/3*y**3)

def Fx (Pot, x, y, h= 1e-3):
    return (Pot(x+h, y) - Pot(x, y)) / h


def Fy (Pot, x, y, h = 1e-3):
    return (Pot(x,y+h) - Pot(x,y)) / h