import numpy as np
np.float = float
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from Potentiel import *
from RK4 import *
from RK2 import *
from Euler import *
import random


def Orbite(wn, N, h, Methode, pot) :
    
    Trajectoire = np.zeros((4,N))

    for i in range(N) : 
        Trajectoire[:,i] = wn
        wn = Methode(wn, f, h,pot)

    return Trajectoire


def Lin_Regression(X,Y) :
    linear_regressor = LinearRegression(fit_intercept=True)
    result = linear_regressor.fit(X, Y)

    return result.coef_[0], result.intercept_


def Gottwald_Melbourne_v1(wn, N, h) :
    
    liste_a = []
    norbite = 0
    for k in range(100):
        Trajectoire = np.zeros((4,N))
        p= np.zeros(N)
        theta= np.zeros(N)
        M= np.zeros(N)
        tau= np.zeros(N)

        #t= np.linspace(0, N*h, N)
        c= 1.7

        Trajectoire[:,0]= wn
        wn = RK4(wn, f, h, Henon_Heiles)

        for ntau in range(1,N) : 

            Trajectoire[:,ntau] = wn
            wn = RK4(wn, f, h, Henon_Heiles)
            tau[ntau]= ntau*h
    

        theta= c*tau + h*np.cumsum(Trajectoire[0]+Trajectoire[1])
        p= h* np.cumsum((Trajectoire[0]+Trajectoire[1])*np.cos(theta))
        #data_M = (p[1:len(p)]-p[0:len(p)-1])**2/(N*h*np.ones(N-1)-tau[0:N-1])
        #M= np.cumsum(data_M)
        
        M = np.zeros(N//10)
        for _ in range(N//10):
            M[_] = 1/(9*N//10) * np.sum((np.roll(p,-_)[0:9*N//10]-p[0:9*N//10])**2)


        a, b = np.polyfit(np.log(tau[0:N//10]+h),np.log(M[0:N//10]+1e-5),1)
        if a < 0.5:
            norbite +=1

    

    #tau= tau.reshape(-1,1)
    #K, b= Lin_Regression(np.log(tau[0:N-4]+1e-7), np.log(M+1))

    #plt.plot(liste_a)
    #plt.show()
    #a = np.mean(np.array([liste_a]))
    """print(a)
    plt.figure()


    plt.scatter(np.log(tau[N//20:N//10]), np.log(M[N//20:N//10]+1e-5))
    plt.show()"""

    #return K

    return norbite /100





def Gottwald_Melbourne_v2(wn, N, h) :
    Trajectoire = np.zeros((4,N))
    p= np.zeros(N)
    #p_mt= np.zeros(N)
    theta= np.zeros(N)
    M= np.zeros(N)
    tau= np.zeros(N)
    K= np.zeros(N)

    #t= np.linspace(0, N*h, N)
    c= 1.7

    Trajectoire[:,0]= wn
    wn = RK4(wn, f, h, Henon_Heiles)

    for ntau in range(1,N) : 

        Trajectoire[:,ntau] = wn
        wn = RK4(wn, f, h, Henon_Heiles)
        print(wn)

        x= Trajectoire[0, 0:ntau+1]
        s= ntau*h

        theta[ntau]= c*s + h*sum(x)
        p[ntau] = h*sum(x*np.cos(theta[0:ntau+1]))
        #p_mt[ntau] = p[ntau-1]

    for ntau in range(1,N) : 
        for nt in range(0,N-ntau) :
            M[ntau]= M[ntau] + (p[nt+ntau]-p[nt])**2
        M[ntau]= M[ntau]/(h*(N-ntau))    
        tau[ntau]= ntau*h
        #print(M[ntau])
    
    tau= tau.reshape(-1,1)
    K, b= Lin_Regression(np.log(tau+1e-7), np.log(M+1))

    """plt.figure()
    plt.scatter(np.log(tau), np.log(M+1))
    plt.show()"""

    return K

def Chaos_Gottwald_Melbourne(N, h, ntraj=300) :
    # energy values for which we will compute the trajectory
    E_values= np.linspace(0.01, 0.16, 14)

    # relative area occupied by the curve in Poincare section
    relative_area= np.zeros(len(E_values))


    for k in range(len(E_values)) : 
        E= E_values[k]
        nb_curve = 0
        e= 0.3

        for n in range(ntraj) :
            b = 0
            # to check that the initial condition will work before computing a trajectory
            while (b==0) : 
                yi = random.uniform(-0.4,0.4)
                vi = random.uniform(-0.4,0.4)

                if (2*(E-Henon_Heiles(0,yi))-vi**2 > 0)  :
                    b = 1

            ui= np.sqrt(2*(E-Henon_Heiles(0,yi))-vi**2)        
            wn= np.array([0,yi,ui,vi])

            K= Gottwald_Melbourne_v1(wn, N, h)
            if K<e :
                nb_curve= nb_curve + 1
                #print(nb_curve)


        # we compute the relative area occupied by the curves 
        relative_area[k]= nb_curve/ntraj

    plt.figure()

    plt.scatter(E_values, relative_area)
    plt.xlabel('Energy')
    plt.ylabel('Relative area')

    plt.show()


    return
            

