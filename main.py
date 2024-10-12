#    _      ______       _____ _    _          ____   _____      _______ ______ _    _ _____         _____   __  _   _ _____          _      
#   | |    |  ____|     / ____| |  | |   /\   / __ \ / ____|  /\|__   __|  ____| |  | |  __ \       / ____|_/_/_| \ | |_   _|   /\   | |     
#   | |    | |__       | |    | |__| |  /  \ | |  | | (___   /  \  | |  | |__  | |  | | |__) |     | |  __| ____|  \| | | |    /  \  | |     
#   | |    |  __|      | |    |  __  | / /\ \| |  | |\___ \ / /\ \ | |  |  __| | |  | |  _  /      | | |_ |  _| | . ` | | |   / /\ \ | |     
#   | |____| |____     | |____| |  | |/ ____ \ |__| |____) / ____ \| |  | |____| |__| | | \ \      | |__| | |___| |\  |_| |_ / ____ \| |____ 
#   |______|______|     \_____|_|  |_/_/    \_\____/|_____/_/    \_\_|  |______|\____/|_|  \_\      \_____|_____|_| \_|_____/_/    \_\______|
#                                                                                                                                            
# 

import numpy as np
import matplotlib.pyplot as plt

from RK2 import*
from Euler import*
from RK4 import*

from Potentiel import*
from Poincare import*
from Energie import*


def Plot_Trajectoires(wn, N, h, pot) :
    """plots the trajectories of a point inside a gravitational potential, computed using Euler, RK2, RK4 methods
    :param wn: vector containing x,y position and u,v velocities
    :param N: number of iterations
    :param h: time step
    :param pot: gravitational potential we consider
    :return: none
    """

    Trajectoire_RK2 = Orbite(wn, N, h, RK2,pot)
    Trajectoire_RK4 = Orbite(wn, N, h, RK4,pot)
    Trajectoire_Euler = Orbite(wn, N, h, Euler,pot)
        
    fig = plt.figure()

    XY = fig.add_subplot(211)
    PXPY = fig.add_subplot(212)

    XY.set_xlabel("X")
    XY.set_ylabel("Y")

    PXPY.set_xlabel("Px")
    PXPY.set_ylabel("Py")

    

    XY.scatter(Trajectoire_Euler[0],Trajectoire_Euler[1],s=1, label="Euler",color="red")
    PXPY.scatter(Trajectoire_Euler[2],Trajectoire_Euler[3],s=1, label="Euler",color="red")


    XY.scatter(Trajectoire_RK2[0],Trajectoire_RK2[1],s=1, label="RK2",color="blue")
    PXPY.scatter(Trajectoire_RK2[2],Trajectoire_RK2[3],s=1, label="RK2",color="blue")

    XY.scatter(Trajectoire_RK4[0],Trajectoire_RK4[1],s=1, label="RK4",color="green")
    PXPY.scatter(Trajectoire_RK4[2],Trajectoire_RK4[3],s=1, label="RK4",color="green")


    XY.legend()
    PXPY.legend()

    plt.show()


   
    


if __name__ == "__main__" :

    # conditions initiales générales
    wn = np.array([0,1,1,0])
    N = int(10 * 10**3)
    h = 10.**-3
    pot = Kepler
    
    # paramètres pour la section de Poincaré (ancienne version)
    
    """error = 1e-3
    E = 1/12
    y= np.arange(-0.4, 0.4, 0.05)           
    v= np.arange(-0.4, 0.4, 0.05)"""   
      
    
    #Plot_Trajectoires(wn, N, h, pot)

    #Poincare_version2(E, y, v, h, error, N, pot)
    
    #Plot_Energie(wn, pot)
    
   

    # -----------------------------------------------------------
    # -----------------------------------------------------------
    # anciens codes :

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
    plt.show()"""
    
    """energies_rk2 = []
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

    wn = np.array([0,1,1,0])

    n = 1

    N = int(10000 * 10**n)
    h = 10.**-n
    pot = Henon_Heiles

    # paramètres supplémentaires pour la section de Poincaré
    E = 1/12
        
    

    #Poincarre_version1(wn, N, h, pot)
    Poincare_version2(E, h, N, pot, ntraj=10)
    
    """energies_rk2 = []
    energies_rk4 = []
    energies_euler = []

    dt = np.arange(1,5,0.5)

    for n in dt : 
        wn = np.array([0,1,1,0])
        N = int(10 * 10**n)
        h = 10.**-n
        pot = Kepler

    if(0) :
        for n in np.arange(1,5,1) : 
            wn = np.array([0,1,1,0])
            N = int(10 * 10**n)
            h = 10.**-n
            pot = Kepler

            #Plot_Trajectoires(wn, N, h, pot)
            #Poincarre (wn, N, h, pot)

            Energie(wn, N, h, pot)
        energies_euler.append(Energie(wn, N, h, pot)[0])
        energies_rk2.append(Energie(wn, N, h, pot)[1])
        energies_rk4.append(Energie(wn, N, h, pot)[2])

    plt.plot(np.log10(10**-dt), energies_euler,label="Euler")
    plt.plot(np.log10(10**-dt), energies_rk2,label="RK2")
    plt.plot(np.log10(10**-dt), energies_rk4,label="RK4")

    axes = plt.gca()

    axes.invert_xaxis()

    axes.set_xlabel(r"$\ln \Delta t [s]$")
    axes.set_ylabel(r"max $\ln \left(\Delta E / E \right)$")

    plt.legend()
    plt.show()
        """
    
    
    """wn = np.array([0,0.1,0,0.097])
    pot = Kepler

    dt_Euler, E_Euler= Energie_std(wn, pot, Euler)
    dt_RK2, E_RK2= Energie_std(wn, pot, RK2)
    dt_RK4, E_RK4= Energie_std(wn, pot, RK4)

    plt.figure()
    plt.plot(dt_Euler, E_Euler, label="Euler")
    plt.plot(dt_RK2, E_RK2, label="RK2")
    plt.plot(dt_RK4, E_RK4, label="RK4")
    plt.ylabel("log(dt), dt in s")
    plt.xlabel("log(stdE/E0), E in J")  # standard deviation
    plt.legend()
    axe= plt.gca()
    axe.invert_xaxis()
    plt.show()"""
    