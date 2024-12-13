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

    # Teste les diffÃ©rentes mÃ©thodes de rÃ©solution de l'Ã©quation du mouvement en comparant la conservation de l'Ã©nergie

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

    liste_poincarres = []

    for _ in range(ntraj):

        liste_poincarres.append(Poincarre_test(E,h,N,pot))

    Poincarre_solver(liste_poincarres,E,h,N,pot,deux=False, plot=True)

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

        liste_poincarres = []

        for j in range(1000):

            liste_poincarres.append(Poincarre_test(E,h,N,pot))

        solver = Poincarre_solver(liste_poincarres,E,h,N,pot,deux=True, plot=False)
        mesure = solver.Chaos_measure(muc=5e-2)
        Resultat_chaos[i] = mesure[0]
        mu_moyen[i] = mesure[1]
        mu_std[i] = mesure[2]

    axes = plt.gca()
    axes.set_xlabel("E [J]")
    axes.set_ylabel("Surface relative")

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

    liste_poincarres = []

    for j in range(1000):

            liste_poincarres.append(Poincarre_test(E_50,h,N,pot))
    solver = Melbourne_solver(liste_poincarres,E_50,h,N,pot, plot=False)
    a_crit = solver.Chaos_measure(calibre = False)
    
    for i in range(len(liste_E)) :
        
        E = liste_E[i]
                
        print("E = "+str(E))

        liste_poincarres = []

        for j in range(1000):

            liste_poincarres.append(Poincarre_test(E,h,N,pot))

        solver = Melbourne_solver(liste_poincarres,E,h,N,pot, plot=False)
        mesure = solver.Chaos_measure(a_crit=a_crit, calibre = True)
        Resultat_chaos[i] = mesure

    axes = plt.gca()
    axes.set_xlabel("E [J]")
    axes.set_ylabel("Surface relative")

    plt.plot(liste_E,Resultat_chaos)
    
    plt.legend()
    plt.show()

    """

    liste_E = np.arange(0.02,0.165,0.02)
    liste_a = []

    for i in range(len(liste_E)) :
        
        E = liste_E[i]
                
        print("E = "+str(E))

        liste_poincarres = []


        p = Poincarre_test(E,h,N,pot)
        wn = np.array([0,p.yi,np.sqrt(2*(E-pot(0,p.yi))-p.vi**2),p.vi])

        liste_a.append(Gottwald_Melbourne_v1(wn,N,h))

    liste_a = np.array(liste_a)
    plt.plot(liste_E, liste_a)

    axes = plt.gca()
    axes.set_xlabel("E [J]")
    #axes.set_ylabel(r"$K = \lim_{+\inf} \log M(t) / \log t$")
    axes.set_ylabel(r"Proportion ergodique")

    plt.show()"""

def test_chaos_1_parall (E, Result_chaos_dict, h=10.**-1, N=1000, pot=Henon_Heiles):

    # Test du chaos Ã  partir de la premiÃ¨re mÃ©thode
    # Affiche une mesure du chaos en fonction de l'Ã©nergie, et la distribution des mus

    yi = []
    vi = []

    Resultat_chaos = 0.
    mu_moyen = 0.
    mu_std = 0.
            
    print("E = "+str(E))

    liste_poincarres = []

    for j in range(1000):

        liste_poincarres.append(Poincarre_test(E,h,N,pot))

    solver = Poincarre_solver(liste_poincarres,E,h,N,pot,deux=True, plot=False)
    mesure = solver.Chaos_measure(muc=5e-2)
    Resultat_chaos = mesure[0]
    mu_moyen = mesure[1]
    mu_std = mesure[2]

    print(Resultat_chaos)
    Result_chaos_dict["result-chaos"]= Resultat_chaos


    return 



if __name__ == "__main__" :

    n = 2

    N = int(4000)
    h = 0.05
    pot = Henon_Heiles
    wn = np.array([0,0.1,0.157,0.1])
    E = 1/16
    
    liste_E = np.arange(0.02,0.165,0.02)
    Result_chaos_dict= {}

    """with Pool() as pool :
        result= pool.map(test_chaos_1_parall, range(0.02, 0.165, 0.02))
        Result_chaos.append(result)

    print("Process finished")"""


    for i in range(len(liste_E)//4) : 
        p1= Process(target=test_chaos_1_parall, args=(liste_E[4*i],Result_chaos_dict,))
        p2= Process(target=test_chaos_1_parall, args=(liste_E[4*i +1],Result_chaos_dict,))
        p3= Process(target=test_chaos_1_parall, args=(liste_E[4*i +2],Result_chaos_dict,))
        p4= Process(target=test_chaos_1_parall, args=(liste_E[4*i +3],Result_chaos_dict,))

        p1.start()
        p2.start()
        p3.start()
        p4.start()

        p1.join()
        p2.join()
        p3.join()
        p4.join()

    Result_chaos= Result_chaos_dict.get("result-chaos")
    print(Result_chaos)
    Energy= np.array(liste_E)

    #plt.figure()
    #plt.scatter(Energy, Result_chaos)
    #plt.show()

    #p1= Process(target=test_chaos_1_parall, args=(0.02,))
    #p2= Process(target=test_chaos_1_parall, args=(0.06,))

    """p1.start()
    print('processing 1')
    p2.start()
    print('processing 2')
    p1.join()
    p2.join()"""


    """liste_E = np.arange(0.02,0.165,0.02)
    table= []

    p= Process(target=test_chaos_1_parall, args=(liste_E))
    p.start()
    p.join()

    for i in range(0, len(liste_E)) :
        table.append([liste_E[i], h, N, "Henon_Heiles"])

    with Pool() as pool :
        result= pool.starmap(test_chaos_1_parall, liste_E)
    print("End program")"""