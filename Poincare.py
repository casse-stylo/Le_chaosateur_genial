import numpy as np
import matplotlib.pyplot as plt

from Potentiel import *
from RK4 import *
import random

def Poincare_version2(E, h, N, Pot, Method=RK4, ntraj = 10) : 

    plt.figure()

    plt.xlabel("y [m]")
    plt.ylabel("v [m/s]")
    plt.title(f"Poincare section for E =  {round(E,3)} J")


    for n in range(ntraj) :
        
        y_list = np.array([])
        v_list = np.array([])

        yi = random.uniform(-0.4,0.4)
        vi = random.uniform(-0.4,0.4)

        if 2*(E-Pot(0,yi))-vi**2 < 0 : continue 

        u= np.sqrt(2*(E-Pot(0,yi))-vi**2)        # initial x velocity
        wn= np.array([0,yi,u,vi])

        # computation of the trajectory and Poincare section
        wnmoins1 = wn   

        for i in range(N) : 
            if np.sign(wnmoins1[0]) != np.sign(wn[0]) :

                # Interpolation linéaire

                y0 = wnmoins1[1] - wnmoins1[0] * (wn[1] - wnmoins1[1]) / (wn[0] - wnmoins1[0])
                v0 = wnmoins1[3] - wnmoins1[0] * (wn[3] - wnmoins1[3]) / (wn[0] - wnmoins1[0])

                y_list = np.append(y_list,y0)
                v_list = np.append(v_list,v0)

            wnmoins1 = wn
            wn = Method(wn, f, h, Pot)
            
        
        #plt.scatter(y_list, v_list, s=0.5, color="red")

    #plt.show()

    return 



