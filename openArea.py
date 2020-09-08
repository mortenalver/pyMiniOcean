from scenario import *
from oceanState import *
import utils
import math
import netcdfStorage

class OpenArea(Scenario):

    os = None
    bounds = {}

    def initialize(self, sp):
        sp.dt = 100#120#240#18#45#15
        sp.nsub = 10
        sp.dx = 4000
        sp.KBi = 1e9

        self.E_floating = False
        self.U_floating = True
        self.V_floating = True

        self.os = sp.getOcean(220, 200, 5)
        self.os.dz[:] = np.array([10, 10, 20, 20, 40])

        self.os.depth[:] = 100
        ## Shallow area:
        self.os.depth[int(0.4*self.os.imax):int(0.6*self.os.imax), \
            int(0.4*self.os.jmax):int(0.6*self.os.jmax)] = self.os.depth[0,0]*0.45
        self.os.depth[...] = utils.smooth(self.os.depth, D=0.1, repeats=5)
        #self.os.depth[:] = netcdfStorage.getDepthMatrix('../matlab/openArea.nc')
        #self.os.depth[:] = netcdfStorage.getDepthMatrix('../matlab/test.nc')
        self.os.calcKmmDzz()

        tval = [12, 11, 10, 9, 8]
        sval = [25, 27, 28, 29, 30]

        for i in range(0,self.os.imax):
            for j in range(0,self.os.jmax):
                #self.os.E[i,j] = 0.2 - 0.4*i/self.os.imax
                for k in range(0,self.os.kmm[i,j]):
                    if i<self.os.imax-1:
                        self.os.U[i,j,k] = -0*0.05
                    self.os.T[i,j,k] = tval[k]
                    self.os.S[i,j,k] = sval[k] #+ 3.5*(i/self.os.imax)*math.exp(-(self.os.midLayerDepths[k]-5)/10)

        self.initBounds(self.os, self.os.imax, self.os.jmax, self.os.kmax, tval, sval)

    def getOs(self):
        return self.os


    def initBounds(self, os, imax, jmax, kmax, tval, sval):

        sleft = np.zeros((jmax,kmax))
        sright = np.zeros((jmax, kmax))
        sbottom = np.zeros((imax, kmax))
        sTop = np.zeros((imax, kmax))
        buleft = np.zeros((jmax, os.kmax))
        buright = np.zeros((jmax, os.kmax))
        for i in range(0, jmax):
            vval = 3*(i/jmax - math.pow(i/jmax, 2))
            if i<jmax/5 or i>4*jmax/5:
                vval = 0.5*vval
            buleft[i,:] = np.ones((kmax))*vval
            buright[i, :] = np.ones((kmax)) * vval
            for k in range(0,kmax):
                sleft[i] = sval[k] + 0.1*vval
                sright[i] = sval[k] + 0.1 * vval
        for i in range(0, imax):
            for k in range(0, kmax):
                sTop[i] = sval[k]
                sbottom[i] = sval[k]

        self.bounds['U'] = [buleft, np.zeros((imax-1, os.kmax)),\
                       buright, np.zeros((imax-1, os.kmax))]
        self.bounds['V'] = [np.zeros((jmax-1, os.kmax)), np.zeros((imax, os.kmax)), \
                      np.zeros((jmax-1, os.kmax)), np.zeros((imax, os.kmax))]
        self.bounds['S'] = [sleft, sbottom, sright, sTop]

    def setBounds(self, imax, jmax ,kmax, fullDepth, t, os):

        return self.bounds