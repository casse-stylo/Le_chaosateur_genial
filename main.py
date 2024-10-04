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

def Orbite(wn, N, Methode = RK2) :
    
    Trajectoire = np.zeros((4,N))

    for i in range(N) : 
        Trajectoire[:,i] = wn
        wn = Methode(wn, f, 1e-3,pot=Kepler)

    return Trajectoire

if __name__ == "__main__" :


    # RK2    
    wn = np.array([1,0,0,1])
    N = 10000

    Trajectoire_RK2 = Orbite(wn, N, RK2)
        
    fig = plt.figure()

    XY = fig.add_subplot(211)
    PXPY = fig.add_subplot(212)

    XY.set_xlabel("X")
    XY.set_ylabel("Y")

    PXPY.set_xlabel("Px")
    PXPY.set_ylabel("Py")

    XY.scatter(Trajectoire_RK2[0],Trajectoire_RK2[1],s=1, label="RK2")
    PXPY.scatter(Trajectoire_RK2[2],Trajectoire_RK2[3],s=1, label="RK2")

    XY.legend()
    PXPY.legend()

    plt.show()