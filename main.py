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

def Orbite(wn, N, h, Methode = RK2) :
    
    Trajectoire = np.zeros((4,N))

    for i in range(N) : 
        Trajectoire[:,i] = wn
        wn = Methode(wn, f, h,pot=Kepler)

    return Trajectoire

if __name__ == "__main__" :


    wn = np.array([1,0,0,1])
    N = 100000
    h = 1e-3

    Trajectoire_RK2 = Orbite(wn, N, h, RK2)
    Trajectoire_RK4 = Orbite(wn, N, h, RK4)

    Trajectoire_Euler = Orbite(wn, N, h, Euler)
        
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