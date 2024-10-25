import numpy as np
import matplotlib.pyplot as plt

from Potentiel import *
from RK4 import *
import random


def Chaos_HenonHeiles(h, N, Pot=Henon_Heiles, muc= 1e-4, Method=RK4, ntraj = 300) : 

    # energy values for which we will compute the trajectory
    E_values= np.linspace(0.01, 0.17, 15)

    # relative area occupied by the curve in Poincare section
    relative_area= np.zeros(len(E_values))


    for k in range(len(E_values)) : 
        E= E_values[k]
        nb_curve = 0

        for n in range(ntraj) :
            # we need to count the number of points in the Poincare section for each trajectory
            count1= 0  
            count2= 0

            mu= 0.
            b = 0
            # to check that the initial condition will work before computing a trajectory
            while (b==0) : 
                yi1 = random.uniform(-0.4,0.4)
                vi1 = random.uniform(-0.4,0.4)

                # second initial point at a distance 1e-8
                yi2= yi1 - np.sqrt(0.5*1e-8)
                vi2= vi1 - np.sqrt(0.5*1e-8)

                if (2*(E-Pot(0,yi1))-vi1**2 > 0) & (2*(E-Pot(0,yi2))-vi2**2 > 0) :
                    b = 1

            
            y1_list = np.array([])
            v1_list = np.array([])

            y2_list = np.array([])
            v2_list = np.array([])

            u1= np.sqrt(2*(E-Pot(0,yi1))-vi1**2)        
            wn1= np.array([0,yi1,u1,vi1])

            u2= np.sqrt(2*(E-Pot(0,yi2))-vi2**2)        
            wn2= np.array([0,yi2,u2,vi2])

            # computation of the trajectory and Poincare section
            wn1_moins1 = wn1 
            wn2_moins1 = wn2


            for i in range(N) : 

                # first trajectory
                while(count1<20) : 
                    if np.sign(wn1_moins1[0]) != np.sign(wn1[0]) :

                        count1= count1+1

                        # Interpolation linéaire

                        y01 = wn1_moins1[1] - wn1_moins1[0] * (wn1[1] - wn1_moins1[1]) / (wn1[0] - wn1_moins1[0])
                        v01 = wn1_moins1[3] - wn1_moins1[0] * (wn1[3] - wn1_moins1[3]) / (wn1[0] - wn1_moins1[0])

                        y1_list = np.append(y1_list,y01)
                        v1_list = np.append(v1_list,v01)

                    wn1_moins1 = wn1
                    wn1 = Method(wn1, f, h, Pot)

                # second trajectory
                while(count2<20) :
                    if np.sign(wn2_moins1[0]) != np.sign(wn2[0]) :

                        count2= count2+1

                        # Interpolation linéaire

                        y02 = wn2_moins1[1] - wn2_moins1[0] * (wn2[1] - wn2_moins1[1]) / (wn2[0] - wn2_moins1[0])
                        v02 = wn2_moins1[3] - wn2_moins1[0] * (wn2[3] - wn2_moins1[3]) / (wn2[0] - wn2_moins1[0])

                        y2_list = np.append(y2_list,y02)
                        v2_list = np.append(v2_list,v02)


                    wn2_moins1 = wn2
                    wn2 = Method(wn2, f, h, Pot)

            # we compute the mu value 
            for j in range(len(y1_list)) :
                mu= mu + (y1_list[j] - y2_list[j])**2 + (v1_list[j] - v2_list[j])**2

            # we check if the trajectory is a curve or ergodic according to mu value
            if mu<muc : 
                nb_curve= nb_curve + 1
                #print(mu)
        
        # we compute the relative area occupied by the curves 
        relative_area[k]= nb_curve/ntraj


    plt.figure()

    plt.scatter(E_values, relative_area)
    plt.xlabel('Energy')
    plt.ylabel('Relative area')

    plt.show()


    return




