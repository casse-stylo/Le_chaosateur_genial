import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*


def RK4 (wn, f, h, pot):
    """apply RK4 method to a 4D phase space vector
    :param wn: vector containing x,y position and u,v velocities
    :param f: function giving the derivatives of the Hamiltonian
    :param h: time step
    :param pot: gravitational potential we consider  
    :result: wn updated according to RK4 scheme
    """
    
    k1 = f(wn, pot)
    k2 = f(wn + h/2 * k1, pot)
    k3 = f(wn+ h/2 * k2, pot)
    k4 = f(wn+h*k3, pot)


    return wn + h/6 * (k1 + 2*k2 + 2*k3 + k4)

