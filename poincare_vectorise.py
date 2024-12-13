import numpy as np
import matplotlib.pyplot as plt
import random

import numpy as np
import matplotlib.pyplot as plt

from Potentiel import *
from RK4 import *
import random
    
    

class Melbourne_solver ():
    def __init__(self, liste_poincarres, E, h, N, Pot, c = 1.7,plot = False,Method=RK4):

        self.E = E
        self.h = h
        self.N = N
        self.Pot = Pot
        

        self.yi = []
        self.vi = []

        self.theta = []
        self.p = []


        for p in liste_poincarres :
            self.yi.append(p.yi)
            self.vi.append(p.vi)

        ntraj = len(liste_poincarres)

        self.liste_poincarre = np.array(liste_poincarres)
        self.epsilon = liste_poincarres[0].epsilon


        self.yi = np.array(self.yi)
        self.vi = np.array(self.vi)

        self.ntraj = len(self.yi)


        u= np.sqrt(2*(self.E-self.Pot(0,self.yi))-self.vi**2)        # initial x velocity

        wn= np.array([[np.zeros(self.ntraj),self.yi,u,self.vi],[np.zeros(self.ntraj),self.yi,u,self.vi]])
        Trajectoires = np.zeros((ntraj,N,4))


        """
        Structure de wn
        premier indice = 0 ou 1 : itÃ©ration prÃ©cÃ©dente ou en cours
        deuxiÃ¨me indice = 0,...,3 : x, y, vx, vy
        troisiÃ¨me indice = 0, ..., ntraj : particule considÃ©rÃ©e
        """

        for _ in range(self.N) : 

            wn[0,:,:] = wn[1,:,:]
            wn[1,:,:] = self.RK4(wn, self.f, self.h, self.Pot)
            Trajectoires[:,_,0] = wn[1,0,:]
            Trajectoires[:,_,1] = wn[1,1,:]            
            Trajectoires[:,_,2] = wn[1,2,:]            
            Trajectoires[:,_,3] = wn[1,3,:]
          
        self.liste_poincarre = liste_poincarres
        self.Trajectoires = Trajectoires

        if plot :
            self.plot()
            
    


    def plot (self, deux = False):

        for p in self.liste_poincarre :
            p.plot(deux)

        plt.show()


    def Chaos_measure(self, a_crit= 0.5, calibre = True) : 

        N = self.N

        nb_curve = 0
        h = self.h

        liste_a = []
        plot = False

        for i in range(len(self.liste_poincarre)) :

        
            #t= np.linspace(0, N*h, N)
            c= 1.7

            Trajectoire = self.Trajectoires[i,:,:]

            tau = np.arange(0,self.N)*self.h


            theta= c*tau + self.h*np.cumsum(Trajectoire[:,0]+Trajectoire[:,1])


            p= self.h* np.cumsum((Trajectoire[:,0]+Trajectoire[:,1])*np.cos(theta))
            #data_M = (p[1:len(p)]-p[0:len(p)-1])**2/(N*h*np.ones(N-1)-tau[0:N-1])
            #M= np.cumsum(data_M)

            M = np.zeros(N//10)
            for _ in range(N//10):
                M[_] = 1/(9*N//10) * np.sum((np.roll(p,-_)[0:9*N//10]-p[0:9*N//10])**2)


            a, b = np.polyfit(np.log(tau[0:N//10]+h),np.log(M[0:N//10]+1e-5),1)
            liste_a.append(a)

            if a < a_crit:
                nb_curve +=1



     
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
        
        liste_a = np.array(liste_a)
        
        
        print(nb_curve)
        print(np.median(liste_a))
        
        if calibre :
            return nb_curve/len(self.liste_poincarre)
        else :
            return np.median(liste_a)
        


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


class Poincarre_solver():
    def __init__(self, liste_poincarres, E, h, N, Pot, plot = False, deux= False, Method=RK4):

        self.E = E
        self.h = h
        self.N = N
        self.Pot = Pot

        self.yi = []
        self.vi = []

        self.yi2 = []
        self.vi2 = []


        for p in liste_poincarres :
            self.yi.append(p.yi)
            self.vi.append(p.vi)

            self.yi2.append(p.yi2)
            self.vi2.append(p.vi2)

        self.liste_poincarre = np.array(liste_poincarres)
        self.epsilon = liste_poincarres[0].epsilon


        self.yi = np.array(self.yi)
        self.vi = np.array(self.vi)

        self.yi2 = np.array(self.yi2)
        self.vi2 = np.array(self.vi2)


        self.ntraj = len(self.yi)


        u= np.sqrt(2*(self.E-self.Pot(0,self.yi))-self.vi**2)        # initial x velocity
        u2 = np.sqrt(2*(self.E-self.Pot(0,self.yi2))-self.vi2**2)

        wn= np.array([[np.zeros(self.ntraj),self.yi,u,self.vi],[np.zeros(self.ntraj),self.yi,u,self.vi]])
        wn2 = np.array([[np.zeros(self.ntraj),self.yi2,u2,self.vi2],[np.zeros(self.ntraj),self.yi2,u2,self.vi2]])

        """
        Structure de wn
        premier indice = 0 ou 1 : itÃ©ration prÃ©cÃ©dente ou en cours
        deuxiÃ¨me indice = 0,...,3 : x, y, vx, vy
        troisiÃ¨me indice = 0, ..., ntraj : particule considÃ©rÃ©e
        """

        for _ in range(self.N) : 

            wn[0,:,:] = wn[1,:,:]
            wn[1,:,:] = self.RK4(wn, self.f, self.h, self.Pot)

            signes = np.sign(wn[1,0,:]) - np.sign(wn[0,0,:]) !=0 

            wn_signe = wn[:,:,signes]

            y0 = wn_signe[0,1,:] - wn_signe[0,0,:] * (wn_signe[1,1,:] - wn_signe[0,1,:]) / (wn_signe[1,0,:] - wn_signe[0,0,:])
            v0 = wn_signe[0,3,:] - wn_signe[0,0,:] * (wn_signe[1,3,:] - wn_signe[0,3,:]) / (wn_signe[1,0,:] - wn_signe[0,0,:])

            indices = np.where(signes)

            if np.shape(y0) !=(0,):

                for s in range(len(y0)):

                    liste_poincarres[indices[0][s]].ylist.append(y0[s])
                    liste_poincarres[indices[0][s]].vlist.append(v0[s])

            if deux :

                wn2[0,:,:] = wn2[1,:,:]
                wn2[1,:,:] = self.RK4(wn2, self.f, self.h, self.Pot)

                signes2 = np.sign(wn2[1,0,:]) - np.sign(wn2[0,0,:]) !=0 

                wn_signe2 = wn2[:,:,signes2]

                y02 = wn_signe2[0,1,:] - wn_signe2[0,0,:] * (wn_signe2[1,1,:] - wn_signe2[0,1,:]) / (wn_signe2[1,0,:] - wn_signe2[0,0,:])
                v02 = wn_signe2[0,3,:] - wn_signe2[0,0,:] * (wn_signe2[1,3,:] - wn_signe2[0,3,:]) / (wn_signe2[1,0,:] - wn_signe2[0,0,:])

                indices2 = np.where(signes2)

                if np.shape(y02) !=(0,):

                    for s in range(len(y02)):
                        #print(liste_poincarres[signes2[s]])

                        liste_poincarres[indices2[0][s]].ylist2.append(y02[s])
                        liste_poincarres[indices2[0][s]].vlist2.append(v02[s])

        self.liste_poincarre = liste_poincarres

        if plot :
            self.plot(deux = deux)


    def plot (self, deux = False):

        for p in self.liste_poincarre :
            p.plot(deux)

        plt.show()

    def Chaos_measure(self, muc= 1e-4) : 

        muc = 130000 * self.epsilon

        nb_curve = 0
        liste_mus = []
        mus = []
        std = []
        
        plot = False

        for p in self.liste_poincarre :
            # we need to count the number of points in the Poincare section for each trajectory

    
            # we compute the mu value 

            index = min([25, len(p.ylist), len(p.ylist2)])

            distances = (np.array(p.ylist[:index])-np.array(p.ylist2[:index]))**2 + (np.array(p.vlist[:index]) - np.array(p.vlist2[:index]))**2


            while len(distances)< 25:
                distances = np.append(distances,0)

            mu = np.sum(distances)
            mus.append(mu)
            #print((np.array(p.ylist[:25])-np.array(p.ylist2[:25]))**2 + (np.array(p.vlist[:25]) - np.array(p.vlist2[:25]))**2)
            #print(mu)

            # we check if the trajectory is a curve or ergodic according to mu value

            if mu<muc : 
                nb_curve= nb_curve + 1

            liste_mus.append(distances)

        if plot :

            liste_mus = np.array(liste_mus)

            #for i in range(100):
            #    plt.scatter(range(1,26),liste_mus[i])

            plt.clf()
            plt.title("E = "+str(self.E))
            plt.plot(np.mean(liste_mus,axis = 0))
            plt.savefig("mu "+str(self.E)+".png")

        # we compute the relative area occupied by the curves 
        
        """plt.scatter(np.arange(len(mus)),mus)
        plt.axhline(np.mean(mus))
        plt.axhline(np.mean(mus)-np.std(mus))
        plt.axhline(np.mean(mus)+np.std(mus))
        plt.show()"""

        return nb_curve/len(self.liste_poincarre), np.mean(mus), np.std(mus)


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
    def __init__ (self, E, h, N, Pot, Method=RK4,epsilon = 1e-8):

        self.E = E
        self.h = h
        self.N = N
        self.Pot = Pot
        self.epsilon = epsilon
        self.Method = Method

        self.ylist = []
        self.vlist = []

        self.ylist2 = []
        self.vlist2 = []

        self.CI()

        #self.poincarre()

        #self.plot()
                
  
    def CI(self):

        self.yi = random.uniform(-0.4, 0.4)
        self.vi = random.uniform(-0.4, 0.4)

        self.yi2= self.yi + np.sqrt(self.epsilon)
        self.vi2= self.vi - np.sqrt(self.epsilon)

        
        if 2*(self.E-self.Pot(0,self.yi))-self.vi**2 < 0 : self.CI()
        if 2*(self.E-self.Pot(0,self.yi2))-self.vi2**2 < 0 : self.CI()

    
    def plot(self, deux = False) :

        self.list=np.array([self.ylist, self.vlist])
        self.list2=np.array([self.ylist2, self.vlist2])

        plt.scatter(self.list[0],self.list[1],s=0.5,color="red")
        plt.xlabel("y")
        plt.ylabel("v")
        if deux : plt.scatter(self.list2[0],self.list2[1],s=0.5,color="blue")

