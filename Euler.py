
import numpy as np
import matplotlib.pyplot as plt

def Kepler(x, y):
    return -1/np.sqrt(x**2 + y**2)

def Fx (Pot, x, y, h= 1e-2):
    return (Pot(x+h, y) - Pot(x, y)) / h


def Fy (Pot, x, y, h = 1e-2):
    return (Pot(x,y+h) - Pot(x,y)) / h


def f(wn, pot=Kepler):
    return np.array([wn[2], wn[3], -Fx(pot,wn[0], wn[1]), -Fy(pot,wn[0], wn[1])])


def Euler(wn, f, h, pot) :
    wn = wn + h*f(wn, pot)
    return wn 

#wn = np.array([1,0,0,1])



if __name__ == "__main__" :
    Nstep=1000
    wn= np.zeros((Nstep, 4))
    wn[0,:]= np.array([1,0,0,1])

    wn= Euler(wn, f, 1e-2, Nstep)

    plt.figure()

    plt.scatter(wn[:,0],wn[:,1])

    plt.show()
