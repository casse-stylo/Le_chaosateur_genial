import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*

#test : cercle

def RK2 (wn, f, h, pot):
    """apply RK2 method to a 4D phase space vector
    :param wn: vector containing x,y position and u,v velocities
    :param f: function giving the derivatives of the Hamiltonian
    :param h: time step
    :param pot: gravitational potential we consider  
    :result: wn updated according to RK2 scheme
    """

    wn12 = wn + h*0.5*f(wn,pot)
    wn12p = f(wn12,pot)

    return wn + h*wn12p

