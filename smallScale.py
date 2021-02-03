from scenario import *
from oceanState import *
import utils
import math

class SmallScale(Scenario):

    os = None
    extent1 = 0.2
    extent2 = 0.6
    uWind = None
    vWind = None
    u_w = 0
    v_w = 0
    dt = 40

    def initialize(self, sp):
        sp.dt = self.dt
        sp.nsub = 20
        sp.dx = 1000
        sp.KBi = 2e7
        self.U_floating = True
        self.V_floating = True

        #self.os = sp.getOcean(40, 20, 5)
        self.os = sp.getOcean(220, 300, 5)
        self.os.dz[:] = np.array([10, 10, 20, 20, 40])
        self.os.depth[:] = 100

        for i in range(0, self.os.imax):
            for j in range(0, self.os.jmax):
                self.os.depth[i,j] = 60+20*(math.pow(math.sin(i*7/self.os.imax),2) +
                                          math.pow(math.sin((j-8)*5 / self.os.jmax), 2))
                #self.os.E[i,j] = self.extent - 2*self.extent*i/self.os.imax


        self.os.calcKmmDzz()



        tval = [12, 11, 10, 9, 8]
        sval = [25, 27, 28, 29, 30]

        for i in range(0,self.os.imax):
            for j in range(0,self.os.jmax):
                for k in range(0,self.os.kmm[i,j]):
                    self.os.T[i,j,k] = tval[k]
                    self.os.S[i,j,k] = sval[k]
                    if i < self.os.imax - 1:
                        self.os.U[i,j,k] = 0.25
                    if j>0.4 * self.os.jmax and j<0.6 * self.os.jmax:
                        self.os.S[i, j, k] = self.os.S[i,j,k] + 2/((k+1)*(k+1))
                        if i<self.os.imax-1:
                            self.os.U[i,j,k] = 1

        self.uWind = np.zeros((self.os.imax-1, self.os.jmax))
        self.vWind = np.zeros((self.os.imax, self.os.jmax-1))

    def getOs(self):
        return self.os


    def setBounds(self, imax, jmax ,kmax, fullDepth, t, os):
        bounds = {}
        eLeft = self.extent1*np.ones(jmax)
        eRight = -self.extent1*np.ones(jmax)
        eBottom = np.ones(imax)
        eTop = np.ones(imax)
        for j in range(int(0.4 * jmax), int(0.6 * jmax)):
            eLeft[j] = self.extent2
            eRight[j] = -self.extent2

        for i in range(0, imax):
            eBottom[i] = self.extent1 - 2*self.extent1*i/imax
            eTop[i] = self.extent1 - 2*self.extent1*i/imax
        bounds['E'] = [eLeft, eBottom, eRight, eTop ]

        uLeft = 0.1*np.ones((jmax,kmax))
        for j in range(int(0.4*jmax),int(0.6*jmax)):
            for k in range(0, kmax):
                uLeft[j,k] = 1
        uRight = uLeft.copy()
        uTop = 0.1 * np.ones((imax-1, kmax))
        uBottom = uTop.copy()

        bounds['U'] = [uLeft, uBottom, uRight, uTop]
        bounds['V'] = [np.zeros((jmax-1, os.kmax)), np.zeros((imax, os.kmax)), \
                      np.zeros((jmax-1, os.kmax)), np.zeros((imax, os.kmax))]
        return bounds

    # The setAtmo method gets dimensions for the full domain as input, and should return
    # matrices containing the U and V components of the wind field.
    def setAtmo(self, imax, jmax, t, os):
        beta = 1/1200
        wind_std = 5
        ff = np.exp(-beta*self.dt)
        self.u_w = ff*self.u_w + np.sqrt(1-ff*ff)*np.random.normal(0, wind_std)
        self.v_w = ff * self.v_w + np.sqrt(1 - ff * ff) * np.random.normal(0, wind_std)
        for i in range(0,imax-1):
            for j in range(0,jmax):
                self.uWind[i,j] = self.u_w
        for i in range(0,imax):
            for j in range(0,jmax-1):
                self.vWind[i,j] = self.v_w

        return self.uWind, self.vWind