from scenario import *
from oceanState import *
import utils
import math

class SmallScale(Scenario):

    os = None
    extent = 0

    def initialize(self, sp):
        sp.dt = 40
        sp.nsub = 20
        sp.dx = 1000
        sp.KBi = 2e7
        self.U_floating = True
        self.V_floating = True

        #self.os = OceanState(60, 50, 5)
        self.os = sp.getOcean(200, 180, 5)
        self.os.dz[:] = np.array([10, 10, 20, 20, 40])
        self.os.depth[:] = 100

        for i in range(0, self.os.imax):
            for j in range(0, self.os.jmax):
                self.os.depth[i,j] = 60+20*(math.pow(math.sin(i*7/self.os.imax),2) +
                                          math.pow(math.sin((j-8)*5 / self.os.jmax), 2))
                self.os.E[i,j] = self.extent - 2*self.extent*i/self.os.imax


        self.os.calcKmmDzz()



        tval = [12, 11, 10, 9, 8]
        sval = [25, 27, 28, 29, 30]

        for i in range(0,self.os.imax):
            for j in range(0,self.os.jmax):
                for k in range(0,self.os.kmm[i,j]):
                    self.os.T[i,j,k] = tval[k]
                    self.os.S[i,j,k] = sval[k]
                    if j<0.5 * self.os.jmax:
                        self.os.S[i, j, k] = self.os.S[i,j,k] + 2/((k+1)*(k+1))


    def getOs(self):
        return self.os


    def setBounds(self, imax, jmax ,kmax, fullDepth, t, os):
        bounds = {}
        eLeft = self.extent*np.ones(jmax)
        eRight = -self.extent*np.ones(jmax)
        eBottom = np.ones(imax)
        eTop = np.ones(imax)
        for i in range(0, imax):
            eBottom[i] = self.extent - 2*self.extent*i/imax
            eTop[i] = self.extent - 2*self.extent*i/imax
        bounds['E'] = [eLeft, eBottom, eRight, eTop ]


        #bounds['U'] = [np.zeros((jmax, os.kmax)), np.zeros((imax-1, os.kmax)),\
        #               np.zeros((jmax, os.kmax)), np.zeros((imax-1, os.kmax))]
        #bounds['V'] = [np.zeros((jmax-1, os.kmax)), np.zeros((imax, os.kmax)), \
        #              np.zeros((jmax-1, os.kmax)), np.zeros((imax, os.kmax))]
        return bounds