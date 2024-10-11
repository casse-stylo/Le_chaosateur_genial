
import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*

def Euler(wn, f, h, pot) :
    """apply Euler method to a 4D phase space vector
    :param wn: vector containing x,y position and u,v velocities
    :param f: function giving the derivatives of the Hamiltonian
    :param h: time step
    :param pot: gravitational potential we consider  
    :result: wn updated according to Euler scheme
    """
    wn = wn + h*f(wn, pot)
    return wn 

