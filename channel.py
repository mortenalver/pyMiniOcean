from scenario import *
from oceanState import *
import utils
import math
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

class Channel(Scenario):

    os = None
    f_par = 0
    dx = 10000
    u_val = 0.25

    def initialize(self, sp):
        sp.dt = 300
        sp.nsub = 20
        sp.dx = self.dx
        sp.KBi = 1e11 # Biharmonic constant

        self.E_floating = True
        #self.U_floating = True
        #self.V_floating = True

        self.os = sp.getOcean(90, 30, 8)
        self.os.dz[:] = np.array([10, 20, 20, 50, 50, 50, 50, 50])
        self.os.depth[:] = 300
        dlevs = [25, 50, 75, 125, 175, 225, 250, 275]
        for i in range(0,self.os.imax):
            lval = 1 #int(1+2*(1.5+math.sin(8*(i-30)/self.os.imax)))
            self.os.depth[i, 0:lval] = 0
            self.os.depth[i, (self.os.jmax-lval):self.os.jmax] = 0
            for j in range(lval, lval+len(dlevs)):
                self.os.depth[i, j] = dlevs[j-lval]
                self.os.depth[i, -(j+1)] = dlevs[j - lval]

        #self.os.depth = utils.smooth(self.os.depth, D=0.075, repeats=5)

        #print(self.os.depth[10,:])
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
                self.os.E[i,j] = 0.2 - ((i-1)/(self.os.imax-1))*0.4
                for k in range(0,self.os.kmm[i,j]):
                    self.os.T[i,j,k] = tval[k]
                    self.os.S[i,j,k] = sval[k]
                    if i < self.os.imax-1:
                        self.os.U[i,j,k] = self.u_val

        # Wind setup:
        self.windU = 0*np.ones((self.os.imax-1, self.os.jmax))
        self.windV = np.zeros((self.os.imax, self.os.jmax-1))


        # Store coriolis parameter for specifying boundaries:
        self.f_par = 2*sp.omega*math.sin(sp.phi0)

    def getOs(self):
        return self.os


    def setBounds(self, imax, jmax ,kmax, fullDepth, t, os):
        bounds = {}


        uLeft =  np.zeros((jmax,kmax))
        uRight = np.zeros((jmax,kmax))
        for j in range(0,jmax):
            if j<=1 or j>=jmax-2:
                uLeft[j] = 0.33*self.u_val
                uRight[j] = 0.33 * self.u_val
            elif j==2 or j==jmax-3:
                uLeft[j] = 0.66 * self.u_val
                uRight[j] = 0.66 * self.u_val
            else:
                uLeft[j] = self.u_val
                uRight[j] = self.u_val

        # # Geostrophic speed - calc dE/dx:
        # de_dx = -self.f_par*self.u_val/9.81
        # e_span = de_dx*jmax*self.dx
        # eLeft = np.zeros(jmax)
        # eRight = np.zeros(jmax)
        # vv = 0.1
        # for j in range(0,jmax):
        #     eLeft[j] = vv*(1-e_span*(0.5 - j/jmax))
        #     eRight[j] = -2*vv + vv*(1-e_span*(0.5-j/jmax))
        #
        # eBottom = np.zeros(imax)
        # eTop = np.zeros(imax)

        #bounds['E'] = [eLeft, eBottom, eRight, eTop ]


        bounds['U'] = [uLeft, np.zeros((imax-1, kmax)),
                        uRight, np.zeros((imax-1, kmax))]
        bounds['V'] = [np.zeros((jmax-1, kmax)), np.zeros((imax, kmax)),
                        np.zeros((jmax-1, kmax)), np.zeros((imax, kmax))]
        return bounds


    def setAtmo(self, imax, jmax, t, os):
        return self.windU, self.windV


