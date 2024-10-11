import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*

from RK2 import*
from Euler import*
from RK4 import*

def Orbite(wn, N, h, Methode, pot) :
    
    Trajectoire = np.zeros((4,N))

    for i in range(N) : 
        Trajectoire[:,i] = wn
        wn = Methode(wn, f, h,pot)

    return Trajectoire

def Energie(wn, N, h, pot, plot = False):

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


# autre fonction qui ne fonctionne pas : 
def Energie_std(wn, pot, Methode) :
    exposant = np.arange(1,6,1)
    dt= 10.**-exposant

    size= np.size(exposant)
    E_std= np.zeros(size)
    

    for n in exposant-1 : 
        h = 10.**-n
        N = int(10 * 10**n)
        E= np.zeros(N)
        for i in range(N) : 
            E[i]= 0.5*(wn[2]**2 + wn[3]**2) + pot(wn[0], wn[1])
            wn = Methode(wn, f, h, pot)
            
        E_std[n]= np.std(E)

    return np.log10(dt), np.log10(E_std/E[0])


def Plot_Energie(wn, pot) :
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