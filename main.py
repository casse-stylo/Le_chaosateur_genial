#    _      ______       _____ _    _          ____   _____      _______ ______ _    _ _____         _____   __  _   _ _____          _      
#   | |    |  ____|     / ____| |  | |   /\   / __ \ / ____|  /\|__   __|  ____| |  | |  __ \       / ____|_/_/_| \ | |_   _|   /\   | |     
#   | |    | |__       | |    | |__| |  /  \ | |  | | (___   /  \  | |  | |__  | |  | | |__) |     | |  __| ____|  \| | | |    /  \  | |     
#   | |    |  __|      | |    |  __  | / /\ \| |  | |\___ \ / /\ \ | |  |  __| | |  | |  _  /      | | |_ |  _| | . ` | | |   / /\ \ | |     
#   | |____| |____     | |____| |  | |/ ____ \ |__| |____) / ____ \| |  | |____| |__| | | \ \      | |__| | |___| |\  |_| |_ / ____ \| |____ 
#   |______|______|     \_____|_|  |_/_/    \_\____/|_____/_/    \_\_|  |______|\____/|_|  \_\      \_____|_____|_| \_|_____/_/    \_\______|
#                                                                                                                                            
# 

import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from Potentiel import*

from RK2 import*
from Euler import*
from RK4 import*
from Poincare import*
from ChaosMeasure import *

import time
from poincare_vectorise import*

def Orbite(wn, N, h, Methode, pot) :
    
    Trajectoire = np.zeros((4,N))

    for i in range(N) : 
        Trajectoire[:,i] = wn
        wn = Methode(wn, f, h,pot)

    return Trajectoire

def Plot_Trajectoires(wn, N, h, pot) :

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

def Poincarre_version1(wn, N, h, pot):
    Trajectoire_RK4 = Orbite(wn, N, h, RK4,pot)

    plt.plot(Trajectoire_RK4[1], Trajectoire_RK4[3])

    axes = plt.gca()
    axes.set_xlabel("Y")
    axes.set_ylabel("V")

    plt.show()


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


if __name__ == "__main__" :

    wn = np.array([0,1,1,0])

    n = 1

    N = int(100 * 10**n)
    h = 0.5*10.**-n
    pot = Henon_Heiles

    # paramètres supplémentaires pour la section de Poincaré
    E = 1/12
        
    
    #Poincarre_version1(wn, N, h, pot)
    #Poincare_version2(E, h, N, pot, ntraj=10)
    

    yi = random.uniform(-0.4, 0.4)
    vi = random.uniform(-0.4, 0.4)

    #Poincarre_version1(wn, N, h, pot)
    #Poincarre_test(E, h, N, pot, ntraj=10)
    
    yi = []
    vi = []

    liste_E = np.arange(0.02,0.165,0.005)
    Resultat_chaos = np.zeros(len(liste_E))
    mu_moyen = np.zeros(len(liste_E))
    mu_std = np.zeros(len(liste_E))

    for i in range(len(liste_E)) :

        print(E)
        
        E = liste_E[i]
        liste_poincarres = []

        for j in range(1000):

            liste_poincarres.append(Poincarre_test(E,h,N,pot))

        solver = Poincarre_solver(liste_poincarres,E,h,N,pot,deux=True, plot=False)
        mesure = solver.Chaos_measure(muc=5e-2)
        Resultat_chaos[i] = mesure[0]
        mu_moyen[i] = mesure[1]
        mu_std = mesure[2]
        print(mesure[1], mesure[2])

    axes = plt.gca()
    axes.set_xlabel("E [J]")
    axes.set_ylabel("Surface relative")

    plt.scatter(liste_E, Resultat_chaos)
    plt.show()

    print(mu_moyen,mu_std)

    plt.errorbar(liste_E,mu_moyen,yerr=mu_std)
    plt.axhline(5e-2,ls="--")
    plt.show()


    # calcul de la chaosité en fonction de l'énergie
    #Chaos_HenonHeiles(liste_poincarres, h, N, ntraj=300)
        
    



    
    """energies_rk2 = []
    energies_rk4 = []
    energies_euler = []

    dt = np.arange(1,5,0.5)

    for n in dt : 
        wn = np.array([0,1,1,0])
        N = int(10 * 10**n)
        h = 10.**-n
        pot = Kepler

        #Plot_Trajectoires(wn, N, h, pot)
        #Poincarre (wn, N, h, pot)

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
    plt.show()"""