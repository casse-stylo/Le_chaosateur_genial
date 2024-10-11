import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*

from RK2 import*
from Euler import*
from RK4 import*


def Energie(wn, N, h, pot, plot = False):
    """computes the total energy along the trajectory for each method
    :param wn: vector containing x,y position and u,v velocities
    :param N: number of iterations
    :param h: time step
    :param pot: gravitational potential we consider 
    :param plot: if True, gives a plot of the energy deviation depending on time for each method
    :return: array of maximum energy deviation for each method
    """

    plt.clf()

    Trajectoire_RK2 = Orbite(wn, N, h, RK2,pot)
    Trajectoire_RK4 = Orbite(wn, N, h, RK4,pot)
    Trajectoire_Euler = Orbite(wn, N, h, Euler,pot)

    Energie_RK2 = (Trajectoire_RK2[2]**2 + Trajectoire_RK2[3]**2)/2 + pot(Trajectoire_RK2[0],Trajectoire_RK2[1])
    Energie_RK4 = (Trajectoire_RK4[2]**2 + Trajectoire_RK4[3]**2)/2 + pot(Trajectoire_RK4[0],Trajectoire_RK4[1])
    Energie_Euler = (Trajectoire_Euler[2]**2 + Trajectoire_Euler[3]**2)/2 + pot(Trajectoire_Euler[0],Trajectoire_Euler[1])

    if plot :
        plt.figure(figsize=(12,8))
        T = np.arange(0,N*h,h)

        plt.plot(T,np.log(np.abs((Energie_Euler - Energie_Euler[0])/Energie_Euler[0])),label="Euler")
        plt.plot(T, np.log(np.abs((Energie_RK2 - Energie_RK2[0])/Energie_RK2[0])), label="RK2")
        plt.plot(T, np.log(np.abs((Energie_RK4 - Energie_RK4[0])/Energie_RK4[0])),label="RK4")

        plt.title(r"$\Delta t = $"+str(h)+" s")

        axes = plt.gca()
        axes.set_xlabel("t [s]")
        axes.set_ylabel(r"$\ln \left( \Delta E / E_0 \right)$")

        plt.legend()

        plt.savefig(str(h)+".png")

    return [np.max(np.log(np.abs((Energie_Euler - Energie_Euler[0])/Energie_Euler[0]))), 
                   np.max(np.log(np.abs((Energie_RK2 - Energie_RK2[0])/Energie_RK2[0]))), 
                   np.max(np.log(np.abs((Energie_RK4 - Energie_RK4[0])/Energie_RK4[0])))]



def Plot_Energie(wn, pot) :
    """plots the maximum energy deviation for each method depending on the time step chosen
    :param wn: vector containing x,y position and u,v velocities
    :param pot: gravitational potential we consider
    :return: none
    """

    energies_rk2 = []
    energies_rk4 = []
    energies_euler = []

    dt = np.arange(1,5,1)
    
    for n in np.arange(1,5,1) : 
            Nn = int(10 * 10**n)
            hn = 10.**-n
            Energie(wn, Nn, hn, pot)
            energies_euler.append(Energie(wn, Nn, hn, pot)[0])
            energies_rk2.append(Energie(wn, Nn, hn, pot)[1])
            energies_rk4.append(Energie(wn, Nn, hn, pot)[2])

    plt.figure()

    plt.plot(dt, energies_euler,label="Euler")
    plt.plot(dt, energies_rk2,label="RK2")
    plt.plot(dt, energies_rk4,label="RK4")

    axes = plt.gca()

    axes.set_xlabel(r"$\ln \Delta t [s]$")
    axes.set_ylabel(r"max $\ln \left(\Delta E / E \right)$")

    plt.legend()
    plt.show()
    return 