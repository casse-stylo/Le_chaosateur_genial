import numpy as np
import matplotlib.pyplot as plt

from Potentiel import *
from RK4 import *

# total energy
E= 1/12 

# potential used
Pot= Henon_Heiles

# initial values
x= 0.                                  # initial position in x
y= np.arange(-0.4, 0.4, 0.05)               # initial position in y
v= np.arange(-0.4, 0.4, 0.05)             # initial y velocity

size_y= np.size(y)
size_v= np.size(v)

plt.figure()

for yi in range(size_y):
    for vi in range(size_v) :
        #print(y[yi]," ", v[vi], "\n")
        if (2*(E-Pot(x,y[yi]))-v[vi]**2 < 0) :
            continue
        else :
            u= np.sqrt(2*(E-Pot(x,y[yi]))-v[vi]**2)        # initial x velocity
            wn= np.array([x,y[yi],u,v[vi]])

            N= 100000   # number of iterations
            h= 1e-2    # time step
            e= 1e-3    # error that we can make on x=0 condition for Poincare section


            # computation of the trajectory and Poincare section

            Trajectoire = np.zeros((4,N))
            Poincare_list= []

            for i in range(N) : 
                if ((wn[0]<e)&(wn[0]>-e))&(wn[2]>0) :
                    Poincare_list.append([wn[1], wn[3]])
                Trajectoire[:,i] = wn
                #print(wn)
                wn = RK4(wn, f, h, Pot)
                


            size= len(Poincare_list)
            Poincare_table= np.zeros((size, 2))
            for i, element in zip(range(size), Poincare_list) :
                Poincare_table[i,:]= np.array(element)


            plt.scatter(Poincare_table[:,0], Poincare_table[:,1], s=0.5, color="red")
            #plt.title(f"y={y[yi]}, v={v[vi]}")

plt.show()

#plt.savefig(fname= r"C:\Documents\Cours université\Astrophysique\M2\Numerical_Methods\fig{E}".format)
