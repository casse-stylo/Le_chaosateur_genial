import numpy as np
import matplotlib.pyplot as plt
from math import*

def Kepler(x, y):
    return -1/sqrt(x**2 + y**2)

def Fx (Pot, x, y, h):
    return (Pot(x+h, y) - Pot(x, y)) / h


def Fy (Pot, x, y, h):
    return (Pot(x,y+h) - Pot(x,y)) / h