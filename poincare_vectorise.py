import numpy as np
import matplotlib.pyplot as plt
import random


import numpy as np
import matplotlib.pyplot as plt

from Potentiel import *
from RK4 import *
import random
    

class Poincarre_solver():
    def __init__(self, liste_poincarres, E, h, N, Pot, Method=RK4):

        self.E = E
        self.h = h
        self.N = N
        self.Pot = Pot

        self.yi = []
        self.vi = []

        for p in liste_poincarres :
            self.yi.append(p.yi)
            self.vi.append(p.vi)

        self.liste_poincarre = np.array(liste_poincarres)

        self.yi = np.array(self.yi)
        self.vi = np.array(self.vi)

        self.ntraj = len(self.yi)


        u= np.sqrt(2*(self.E-self.Pot(0,self.yi))-self.vi**2)        # initial x velocity

        wn= np.array([[np.zeros(self.ntraj),self.yi,u,self.vi],[np.zeros(self.ntraj),self.yi,u,self.vi]])


        """
        Structure de wn
        premier indice = 0 ou 1 : itération précédente ou en cours
        deuxième indice = 0,...,3 : x, y, vx, vy
        troisième indice = 0, ..., ntraj : particule considérée
        """

        for i in range(self.N) : 

            wn[0,:,:] = wn[1,:,:]
            wn[1,:,:] = self.RK4(wn, self.f, self.h, self.Pot)

            signes = np.sign(wn[1,0,:]) - np.sign(wn[0,0,:]) !=0 

            wn_signe = wn[:,:,signes]

            y0 = wn_signe[0,1,:] - wn_signe[0,0,:] * (wn_signe[1,1,:] - wn_signe[0,1,:]) / (wn_signe[1,0,:] - wn_signe[0,0,:])
            v0 = wn_signe[0,3,:] - wn_signe[0,0,:] * (wn_signe[1,3,:] - wn_signe[0,3,:]) / (wn_signe[1,0,:] - wn_signe[0,0,:])

            if np.shape(y0) !=(0,):

                for s in range(len(y0)):

                    liste_poincarres[signes[s]].ylist.append(y0[s])
                    liste_poincarres[signes[s]].vlist.append(v0[s])

        for p in liste_poincarres :
            p.plot()

            """           
            y0 = wnmoins1[1] - wnmoins1[0] * (wn[1] - wnmoins1[1]) / (wn[0] - wnmoins1[0])
            v0 = wnmoins1[3] - wnmoins1[0] * (wn[3] - wnmoins1[3]) / (wn[0] - wnmoins1[0])
            """

            #plt.scatter(y0,v0)

        plt.show()

    def RK4 (self,wn, f, h, pot):
        """apply RK4 method to a 4D phase space vector
        :param wn: vector containing x,y position and u,v velocities
        :param f: function giving the derivatives of the Hamiltonian
        :param h: time step
        :param pot: gravitational potential we consider  
        :result: wn updated according to RK4 scheme
        """
        
        k1 = f(wn[1,:,:], pot)
        k2 = f(wn[1,:,:] + h/2 * k1, pot)
        k3 = f(wn[1,:,:]+ h/2 * k2, pot)
        k4 = f(wn[1,:,:]+h*k3, pot)


        return wn[1,:,:] + h/6 * (k1 + 2*k2 + 2*k3 + k4)


    def Fx (self,Pot, x, y, h= 1e-3):
        """derivation with 5 points
        :param Pot: function we want to derive 
        :param x: position along x axis
        :param y: position along y axis
        :param h: time step
        :return: partial derivative of the function with respect to x
        """
        return 1/(12*h) * (-Pot(x+2*h, y)+ 8*Pot(x+h, y) - 8*Pot(x-h, y) + Pot(x-2*h, y))


    def Fy (self,Pot, x, y, h = 1e-3):
        """derivation with 5 points
        :param Pot: function we want to derive 
        :param x: position along x axis
        :param y: position along y axis
        :param h: time step
        :return: partial derivative of the function with respect to y
        """
        return 1/(12*h) * (-Pot(x,y+2*h)+ 8*Pot(x,y+h) - 8*Pot(x,y-h) + Pot(x,y-2*h))


    # Differential equation function

    def f(self,wn, pot=Kepler):
        """function giving the derivatives of the Hamiltonian
        :param wn: vector containing x,y position and u,v velocities
        :param pot: gravitational potential we consider
        :return: array of Hamiltonian's derivatives 
        """
        return np.array([wn[2], wn[3], -self.Fx(pot,wn[0], wn[1]), -self.Fy(pot,wn[0], wn[1])])


class Poincarre_test ():
    def __init__ (self, E, h, N, Pot, Method=RK4):

        self.E = E
        self.h = h
        self.N = N
        self.Pot = Pot
        self.Method = Method

        self.ylist = []
        self.vlist = []

        self.CI()

        #self.poincarre()

        #self.plot()
                
  
    def CI(self):

        self.yi = random.uniform(-0.4, 0.4)
        self.vi = random.uniform(-0.4, 0.4)
        
        if 2*(self.E-self.Pot(0,self.yi))-self.vi**2 < 0 : self.CI()


    
    def plot(self) :

        self.list=np.array([self.ylist, self.vlist])
        plt.scatter(self.list[0],self.list[1],s=0.5,color="red")
   