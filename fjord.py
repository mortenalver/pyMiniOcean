from scenario import *
from oceanState import *
from oceanStateSplit import *
import utils
import math
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import random

class Fjord(Scenario):

    os = None

    def initialize(self, sp):
        sp.dt = 100
        sp.nsub = 20
        sp.dx = 4000
        sp.KBi = 1e10#1e11 # Biharmonic constant

        self.U_floating = False
        self.V_floating = False

        imax = 130
        jmax = 100
        self.os = sp.getOcean(imax, jmax, 8)
        self.os.dz[:] = np.array([5, 5, 10, 10, 20, 50, 50, 50])
        center = (int(imax/1.8), int(jmax/2))

        self.atmoUpdateInt = 3
        self.atmoUpdateCount = -1
        self.atmoUpdateFactor = 0.15
        self.atmoVec = [5, 0]

        maxdepth = 200
        entrDepth = 80
        entrCoord = (int(jmax*0.23), int(jmax*0.77))
        distf = jmax/2.2


        self.randomgen = random.Random()
        self.randomgen.seed(100)

        shallower = ([150, 25, imax*0.6, jmax*0.65],
                     [110, 40, imax*0.45, jmax*0.56],
                     [100, 34, imax *0.8, jmax * 0.53],
                     [75, 20, imax*0.5, jmax*0.3],
                     [50, 20, imax * 0.6, jmax * 0.32],
                     [50, 20, imax * 0.7, jmax * 0.27])

        dlevs = [35, 70, 115, 150, 200, 250, 300, 350, 375]
        for i in range(0,imax):
            for j in range(0,jmax):
                dist = math.sqrt((i-center[0])*(i-center[0]) + (j-center[1])*(j-center[1]))
                self.os.depth[i,j] = max(0, maxdepth*(1 - math.pow(dist/distf, 4)))
                if self.os.depth[i,j] > 0:
                    for shl in shallower:
                        dist2 = math.sqrt((i-shl[2])*(i-shl[2]) + (j-shl[3])*(j-shl[3]))
                        efac = shl[1]/imax
                        self.os.depth[i,j] = self.os.depth[i,j] - shl[0]*(1-math.exp(efac*dist2)/(1+math.exp(efac*dist2)))

                #if self.os.depth[i,j] > 20:
                #    self.os.depth[i,j] = max(0, self.os.depth[i,j] + self.randomgen.gauss(0, 5))
                if i<center[0] and j>=entrCoord[0] and j<entrCoord[1]:
                    dFromLand = min(j-entrCoord[0], entrCoord[1]-j)
                    entrDepthHere = entrDepth*dFromLand/(dFromLand + (entrCoord[1]-entrCoord[0])/6)
                    self.os.depth[i,j] = max(self.os.depth[i,j], entrDepthHere)

        #plt.figure()
        #plt.pcolor(np.transpose(self.os.depth)), plt.colorbar()
        #plt.show()

        self.os.calcKmmDzz()

        profDepths = [0, 10, 20, 30, 40, 100, 200, 1000]
        #temp = [13, 12, 11, 10, 9, 8, 7, 6]
        #salt = [20, 25, 30, 31, 32, 33, 34, 34.5]
        temp = [10, 10, 10, 10, 10, 10, 10, 10]
        salt = [30, 30, 30, 30, 30, 30, 30, 30]

        f = interp1d(profDepths, temp)
        tval = f(self.os.midLayerDepths)
        f = interp1d(profDepths, salt)
        sval = f(self.os.midLayerDepths)
        print(tval)
        print(sval)

        for i in range(0,self.os.imax):
            for j in range(0,self.os.jmax):
                self.os.E[i,j] = 0*(-0.1+0.2*i/self.os.imax)
                for k in range(0,self.os.kmm[i,j]):
                    self.os.T[i,j,k] = tval[k]
                    self.os.S[i,j,k] = sval[k]


        # River setup:
        j = int(self.os.jmax/2)
        i = imax-1
        while self.os.depth[i,j] == 0:
            i=i-1
        self.rpos = (i, j)
        print(self.rpos)

    def getOs(self):
        return self.os

    #def initPassiveTracer(self):
    #    poss = (int(self.os.imax/2), int(self.os.jmax/2))
    #    self.os.X[poss[0], poss[1],0] = 1
    #    return

    def setBounds(self, imax, jmax ,kmax, fullDepth, t, os):
        bounds = {}
        e_val = 2*0.66*math.sin(2*3.14*t/3600/12.2) + 0.33*math.sin(2*3.14*(t-1800)/3600/12)
        u_val = 0.25*(0.66*math.cos(2*3.14*t/3600/12.2) + 0.33*math.cos(2*3.14*(t-1800)/3600/12))

        eLeft = e_val*np.ones(jmax)
        eRight = 0*np.ones(jmax)
        eBottom = 0*np.ones(imax)
        eTop = 0*np.ones(imax)
        bounds['E'] = [eLeft, eBottom, eRight, eTop ]

        #uLeft = u_val*np.ones((jmax,kmax))
        #uRight = 0*np.ones((jmax,kmax))
        #uBottom = 0*np.ones((imax-1,kmax))
        #uTop = 0*np.ones((imax-1,kmax))
        #bounds['U'] = [uLeft, uBottom, uRight, uTop ]

        vLeft = 0*np.ones((jmax-1,kmax))
        vRight = 0*np.ones((jmax-1,kmax))
        vBottom = 0*np.ones((imax,kmax))
        vTop = 0*np.ones((imax,kmax))
        bounds['V'] = [vLeft, vBottom, vRight, vTop ]


        xLeft = np.ones((jmax,kmax))
        xRight = 0*np.ones((jmax,kmax))
        xBottom = 0*np.ones((imax,kmax))
        xTop = 0*np.ones((imax,kmax))
        bounds['X'] = [xLeft, xBottom, xRight, xTop]

        #bounds['U'] = [np.zeros((fullDims[1], os.kmax)), np.zeros((fullDims[0]-1, os.kmax)),\
        #               np.zeros((fullDims[1], os.kmax)), np.zeros((fullDims[0]-1, os.kmax))]
        #bounds['V'] = [np.zeros((fullDims[1]-1, os.kmax)), np.zeros((fullDims[0], os.kmax)), \
        #               np.zeros((fullDims[1]-1, os.kmax)), np.zeros((fullDims[0], os.kmax))]
        return bounds


    def setAtmo(self, imax, jmax, t, os):
        if self.atmoUpdateCount < 0:
            random.seed(100)
            self.windU = self.atmoVec[0] * np.ones((imax - 1, jmax))
            self.windV = self.atmoVec[1] * np.ones((imax, jmax - 1))
        self.atmoUpdateCount = self.atmoUpdateCount + 1
        if self.atmoUpdateCount == self.atmoUpdateInt:
            self.atmoUpdateCount = 0
            self.atmoVec[0] = (1-self.atmoUpdateFactor)*self.atmoVec[0]+self.atmoUpdateFactor*self.randomgen.gauss(0,13)
            self.atmoVec[1] = (1-self.atmoUpdateFactor)*self.atmoVec[1]+self.atmoUpdateFactor*self.randomgen.gauss(0,13)
            self.windU = self.atmoVec[0] * np.ones((imax - 1, jmax))
            self.windV = self.atmoVec[1] * np.ones((imax, jmax - 1))


        return self.windU, self.windV

    def getFreshwater(self, imax, jmax, t, os):
        return [self.rpos,], ["east",], [(2000, 5, 0),]


