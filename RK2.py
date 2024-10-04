import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*

#test : cercle

def RK2 (wn, f, h, pot):

    wn12 = wn + h*0.5*f(wn,pot)
    wn12p = f(wn12,pot)

    return wn + h*wn12p

