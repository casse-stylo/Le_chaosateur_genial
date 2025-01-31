import numpy as np
np.float = float

import matplotlib.pyplot as plt


# Resolution of the differential equation

def Orbite(wn, N, h, Methode, pot) :
    """computes the trajectory of a point in a gravitational potential over time
    :param wn: vector containing initial x,y position and u,v velocities
    :param N: number of iterations
    :param h: time step
    :param Methode: method used to compute the trajectory
    :param pot: gravitational potential we consider
    :return: array containing the trajectory points in phase space
    """
    
    Trajectoire = np.zeros((4,N))

    for i in range(N) : 
        Trajectoire[:,i] = wn
        wn = Methode(wn, f, h,pot)

    return Trajectoire


# Potentials

def Kepler(x, y):
    """Kepler potential function
    :param x: position along x axis
    :param y: position along y axis
    :return: Kepler potential in (x,y)
    """
    return -1/np.sqrt(x**2 + y**2)

def Henon_Heiles (x,y):
    """Henon-Heiles potential function
    :param x: position along x axis
    :param y: position along y axis
    :return: Henon-Heiles potential in (x,y)
    """
    return 0.5*(x**2 + y**2+2*x**2*y-2/3*y**3)


# Derivatives

def Fx (Pot, x, y, h= 1e-3):
    """derivation with 5 points
    :param Pot: function we want to derive 
    :param x: position along x axis
    :param y: position along y axis
    :param h: time step
    :return: partial derivative of the function with respect to x
    """
    return 1/(12*h) * (-Pot(x+2*h, y)+ 8*Pot(x+h, y) - 8*Pot(x-h, y) + Pot(x-2*h, y))


def Fy (Pot, x, y, h = 1e-3):
    """derivation with 5 points
    :param Pot: function we want to derive 
    :param x: position along x axis
    :param y: position along y axis
    :param h: time step
    :return: partial derivative of the function with respect to y
    """
    return 1/(12*h) * (-Pot(x,y+2*h)+ 8*Pot(x,y+h) - 8*Pot(x,y-h) + Pot(x,y-2*h))


# Differential equation function

def f(wn, pot=Kepler):
    """function giving the derivatives of the Hamiltonian
    :param wn: vector containing x,y position and u,v velocities
    :param pot: gravitational potential we consider
    :return: array of Hamiltonian's derivatives 
    """
    return np.array([wn[2], wn[3], -Fx(pot,wn[0], wn[1]), -Fy(pot,wn[0], wn[1])])



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

