#    _      ______       _____ _    _          ____   _____      _______ ______ _    _ _____         _____   __  _   _ _____          _      
#   | |    |  ____|     / ____| |  | |   /\   / __ \ / ____|  /\|__   __|  ____| |  | |  __ \       / ____|_/_/_| \ | |_   _|   /\   | |     
#   | |    | |__       | |    | |__| |  /  \ | |  | | (___   /  \  | |  | |__  | |  | | |__) |     | |  __| ____|  \| | | |    /  \  | |     
#   | |    |  __|      | |    |  __  | / /\ \| |  | |\___ \ / /\ \ | |  |  __| | |  | |  _  /      | | |_ |  _| | . ` | | |   / /\ \ | |     
#   | |____| |____     | |____| |  | |/ ____ \ |__| |____) / ____ \| |  | |____| |__| | | \ \      | |__| | |___| |\  |_| |_ / ____ \| |____ 
#   |______|______|     \_____|_|  |_/_/    \_\____/|_____/_/    \_\_|  |______|\____/|_|  \_\      \_____|_____|_| \_|_____/_/    \_\______|
#                                                                                                                                            
# 

import numpy as np
np.float = float

import matplotlib.pyplot as plt
from Potentiel import*

from RK2 import*
from Euler import*
from RK4 import*
from Poincare import*
from ChaosMeasure import *

import time
from poincare_vectorise import*
from multiprocessing import Pool
from multiprocessing import Process

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


def test_solvers ():

    # Tests the different integrators for the resolution of the equation of motion in the case of the Keplerian potential, and compare the energy conservation

    energies_rk2 = []
    energies_rk4 = []
    energies_euler = []

    dt = np.arange(1,5,0.5)

    for n in dt : 
        wn = np.array([0,1,1,0])
        N = int(10 * 10**n)
        h = 10.**-n
        pot = Kepler

        energies_euler.append(Energie(wn, N, h, pot)[0])
        energies_rk2.append(Energie(wn, N, h, pot)[1])
        energies_rk4.append(Energie(wn, N, h, pot)[2])

    plt.plot(np.log10(10.**-dt), energies_euler,label="Euler")
    plt.plot(np.log10(10.**-dt), energies_rk2,label="RK2")
    plt.plot(np.log10(10.**-dt), energies_rk4,label="RK4")

    axes = plt.gca()

    axes.invert_xaxis()

    axes.set_xlabel(r"$\ln \Delta t [s]$")
    axes.set_ylabel(r"max $\ln \left(\Delta E / E \right)$")

    plt.legend()
    plt.show()

def test_section_poincarre(E, h, N, pot, ntraj = 50):

    # Trace les sections de poincarre

    liste_particules = []

    for _ in range(ntraj):

        liste_particules.append(Particule(E,h,N,pot))

    Poincarre_solver(liste_particules,E,h,N,pot,deux=False, plot=True)

def test_chaos_1 (h, N, pot):

    # Test du chaos Ã  partir de la premiÃ¨re mÃ©thode
    # Affiche une mesure du chaos en fonction de l'Ã©nergie, et la distribution des mus

    yi = []
    vi = []

    liste_E = np.arange(0.02,0.165,0.02)
    Resultat_chaos = np.zeros(len(liste_E))
    mu_moyen = np.zeros(len(liste_E))
    mu_std = np.zeros(len(liste_E))

    for i in range(len(liste_E)) :
        
        E = liste_E[i]
                
        print("E = "+str(E))

        liste_particules = []

        for j in range(1000):

            liste_particules.append(Particule(E,h,N,pot))

        solver = Poincarre_solver(liste_particules,E,h,N,pot,deux=True, plot=False)
        mesure = solver.Chaos_measure(muc=5e-2)
        Resultat_chaos[i] = mesure[0]
        mu_moyen[i] = mesure[1]
        mu_std[i] = mesure[2]

    axes = plt.gca()
    axes.set_xlabel("Energy")
    axes.set_ylabel("Relative area")

    plt.scatter(liste_E, Resultat_chaos)
    plt.show()

    axes = plt.gca()
    axes.set_xlabel("E [J]")
    axes.set_ylabel(r"$\left<\mu\right> \pm \sigma_{\mu}$")

    plt.errorbar(liste_E,mu_moyen,yerr=mu_std)
    plt.axhline(5e-2,ls="--",label = r"$\mu_c$")

    plt.legend()
    plt.show()

def test_Melbourne (N,h,pot=Henon_Heiles, E_50 = 0.14):

    yi = []
    vi = []
    
    # E_50 permet de calibrer a_crit pour que la moitie des orbites de cette energie soient considerees chaotiques

    liste_E = np.arange(0.02,0.165,0.02)
    Resultat_chaos = np.zeros(len(liste_E))
    mu_moyen = np.zeros(len(liste_E))
    mu_std = np.zeros(len(liste_E))

    liste_particules = []

    for j in range(100):

            liste_particules.append(Particule(E_50,h,N,pot))
    solver = Melbourne_solver(liste_particules,E_50,h,N,pot, plot=False)
    a_crit = solver.Chaos_measure(calibre = False)
    
    for i in range(len(liste_E)) :
        
        E = liste_E[i]
                
        print("E = "+str(E))

        liste_particules = []

        for j in range(100):

            liste_particules.append(Particule(E,h,N,pot))

        solver = Melbourne_solver(liste_particules,E,h,N,pot, plot=False)
        mesure = solver.Chaos_measure(a_crit=a_crit, calibre = True)
        Resultat_chaos[i] = mesure

    axes = plt.gca()
    axes.set_xlabel("E [J]")
    axes.set_ylabel("Surface relative")

    plt.plot(liste_E,Resultat_chaos)
    
    plt.legend()
    plt.show()



if __name__ == "__main__" :

    n = 2

    N = int(8000)
    h = 0.05
    pot = Henon_Heiles
    wn = np.array([0,0.1,0.157,0.1])
    E = 1/16

    #Plot_Trajectoires(wn, N, h, Kepler)
    #test_solvers()
    #test_section_poincarre(E = 1/100, h = 1e-3, N = 100000, pot = Henon_Heiles)    # Choisir N assez grand pour avoir des orbites fermÃ©es
    
    
    #test_chaos_1(h = 1e-1, N = 1000, pot = Henon_Heiles)                             # Choisir N plus petit pour avoir 25 points dans la section de poincarrÃ©
    
    test_Melbourne(N, h)
    #print(Gottwald_Melbourne_v1(wn, N, h))
    #Chaos_Gottwald_Melbourne(N, h, 50)
    



