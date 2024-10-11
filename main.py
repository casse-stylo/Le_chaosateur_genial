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


def Energie(wn, N, h, pot):

    plt.clf()

    plt.figure(figsize=(12,8))

    Trajectoire_RK2 = Orbite(wn, N, h, RK2,pot)
    Trajectoire_RK4 = Orbite(wn, N, h, RK4,pot)
    Trajectoire_Euler = Orbite(wn, N, h, Euler,pot)

    Energie_RK2 = (Trajectoire_RK2[2]**2 + Trajectoire_RK2[3]**2)/2 + pot(Trajectoire_RK2[0],Trajectoire_RK2[1])
    Energie_RK4 = (Trajectoire_RK4[2]**2 + Trajectoire_RK4[3]**2)/2 + pot(Trajectoire_RK4[0],Trajectoire_RK4[1])
    Energie_Euler = (Trajectoire_Euler[2]**2 + Trajectoire_Euler[3]**2)/2 + pot(Trajectoire_Euler[0],Trajectoire_Euler[1])

    
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


if __name__ == "__main__" :


    wn = np.array([0,0.1,0,0.097])
    N = 100000
    h = 1e-3
    pot = Henon_Heiles

    # paramètres supplémentaires pour la section de Poincaré
    error = 1e-3
    E = 1/12
    y= np.arange(-0.4, 0.4, 0.05)           
    v= np.arange(-0.4, 0.4, 0.05)             
    

    #Plot_Trajectoires(wn, N, h, pot)
    #Poincarre_version1(wn, N, h, pot)
    Poincare_version2(E, y, v, h, error, N, pot)

    if(0) :
        for n in np.arange(1,5,1) : 
            wn = np.array([0,1,1,0])
            N = int(10 * 10**n)
            h = 10.**-n
            pot = Kepler

            #Plot_Trajectoires(wn, N, h, pot)
            #Poincarre (wn, N, h, pot)

            Energie(wn, N, h, pot)

    