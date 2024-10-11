import numpy as np
import matplotlib.pyplot as plt

from Potentiel import *
from RK4 import *

def Poincare_version2(E, y, v, h, error, N, Pot, Method=RK4) : 
    x= 0.
    size_y= np.size(y)
    size_v= np.size(v)

    plt.figure()

    plt.xlabel("y")
    plt.ylabel("v")
    plt.title(f"Poincare section for energy {E}")

    for yi in range(size_y):
        for vi in range(size_v) :
            #print(y[yi]," ", v[vi], "\n")
            if (2*(E-Pot(x,y[yi]))-v[vi]**2 < 0) :
                continue
            else :
                u= np.sqrt(2*(E-Pot(x,y[yi]))-v[vi]**2)        # initial x velocity
                wn= np.array([x,y[yi],u,v[vi]])

                # computation of the trajectory and Poincare section
                Trajectoire = np.zeros((4,N))
                Poincare_list= []

                for i in range(N) : 
                    if ((wn[0]<error)&(wn[0]>-error))&(wn[2]>0) :
                        Poincare_list.append([wn[1], wn[3]])
                    Trajectoire[:,i] = wn
                    wn = Method(wn, f, h, Pot)
                    


                size= len(Poincare_list)
                Poincare_table= np.zeros((size, 2))
                for i, element in zip(range(size), Poincare_list) :
                    Poincare_table[i,:]= np.array(element)


                plt.scatter(Poincare_table[:,0], Poincare_table[:,1], s=0.5, color="red")

    plt.show()
    return 



