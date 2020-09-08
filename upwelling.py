from scenario import *
from oceanState import *
import utils
import math
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class Upwelling(Scenario):

    river = True
    os = None
    f_par = 0
    dx = 10000
    u_val = -0.5

    def initialize(self, sp):
        sp.dt = 300
        sp.nsub = 20
        sp.dx = self.dx
        sp.KBi = 1e11 # Biharmonic constant

        #self.U_floating = True
        #self.V_floating = True

        self.os = sp.getOcean(60, 52, 10)
        self.os.dz[:] = np.array([5, 5, 10, 10, 20, 50, 50, 50, 100, 200])
        self.os.depth[:] = 400
        dlevs = [35, 70, 115, 150, 200, 250, 300, 350, 375]
        for i in range(0,self.os.imax):
            lval = int(1+2*(1.5+math.sin(8*(i-30)/self.os.imax)))
            self.os.depth[i, 0:lval] = 0
            for j in range(lval, lval+len(dlevs)):
                self.os.depth[i, j] = dlevs[j-lval]
        #self.os.depth[:, -1] = 0
        #self.os.depth[:, -2] = 35
        #self.os.depth[:, -3] = 70
        #self.os.depth[:, -4] = 115
        ##self.os.depth[:, -5] = 150
        #self.os.depth[:, -6] = 200
        #self.os.depth[:, -7] = 250
        #self.os.depth[:, -8] = 300
        #self.os.depth[:, -9] = 350
        #self.os.depth[:, -10] = 375
        self.os.depth = utils.smooth(self.os.depth, D=0.075, repeats=5)

        #plt.figure()
        #plt.pcolor(np.transpose(self.os.depth)), plt.colorbar()
        #plt.show()

        self.os.calcKmmDzz()

        profDepths = [0, 10, 20, 30, 40, 100, 200, 1000]
        temp = [13, 12, 11, 10, 9, 8, 7, 6]
        salt = [20, 25, 30, 31, 32, 33, 34, 34.5]

        f = interp1d(profDepths, temp)
        tval = f(self.os.midLayerDepths)
        f = interp1d(profDepths, salt)
        sval = f(self.os.midLayerDepths)
        print(tval)
        print(sval)

        for i in range(0,self.os.imax):
            for j in range(0,self.os.jmax):
                self.os.E[i,j] = 0
                for k in range(0,self.os.kmm[i,j]):
                    self.os.T[i,j,k] = tval[k]
                    self.os.S[i,j,k] = sval[k]
                    if i < self.os.imax-1:
                        self.os.U[i,j,k] = self.u_val

        # Wind setup:
        self.windU = -10*np.ones((self.os.imax-1, self.os.jmax))
        self.windV = np.zeros((self.os.imax, self.os.jmax-1))

        # River setup:
        if self.river:
            i = int(self.os.imax/2)
            j=0
            while self.os.depth[i,j] == 0:
                j = j+1
            self.rpos = (i, j)
            print(self.rpos)

        # Store coriolis parameter for specifying boundaries:
        self.f_par = 2*sp.omega*math.sin(sp.phi0)

    def getOs(self):
        return self.os


    def setBounds(self, imax, jmax ,kmax, fullDepth, t, os):
        bounds = {}


        uLeft =  self.u_val*np.ones((jmax,kmax))
        uRight = self.u_val*np.ones((jmax,kmax))
        # Geostrophic speed - calc dE/dx:
        de_dx = -self.f_par*self.u_val/9.81
        e_span = de_dx*jmax*self.dx
        e_val = 0*0.1
        eLeft = np.zeros(jmax)
        eRight = np.zeros(jmax)
        for j in range(0,jmax):
            eLeft[j] = -e_span*(0.5 - j/jmax)
            eRight[j] = -e_span*(0.5-j/jmax)

        eBottom = np.zeros(imax)
        eTop = 0.5*e_span*np.ones(imax)
        #for i in range(0, imax):
        #    eBottom[i] = -e_val*(1- 2*i/imax)
        #    eTop[i] = -e_val*(1- 2*i/imax)
        bounds['E'] = [eLeft, eBottom, eRight, eTop ]


        bounds['U'] = [uLeft, np.zeros((imax-1, kmax)),
                       uRight, np.zeros((imax-1, kmax))]
        bounds['V'] = [np.zeros((jmax-1, kmax)), np.zeros((imax, kmax)),
                       np.zeros((jmax-1, kmax)), np.zeros((imax, kmax))]
        return bounds


    def setAtmo(self, imax, jmax, t, os):
        return self.windU, self.windV

    def getFreshwater(self, imax, jmax, t, os):
        if self.river:
            return [self.rpos,], ["south",], [(1000, 8, 0),]
        else:
            return [], [], []

