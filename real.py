from scenario import *
from oceanState import *
import utils
import math
import netcdfStorage

class Real(Scenario):

    os = None

    def initialize(self, sp):
        sp.dt = 1200
        sp.nsub = 25
        sp.dx = 20000
        sp.KBi = 1.5e12  # Biharmonic constant

        self.E_floating = False
        self.U_floating = True
        self.V_floating = True

        #self.os = netcdfStorage.loadState(sp, 'data/real_init.nc', 0)
        self.os = netcdfStorage.loadSINMODState(sp, 'C:/SINMOD_operative/gin/gin_PhysStates.nc', 0, [25, 85, 25, 80])
        #self.os = OceanState(300, 235, 5)
        #self.os.dz[:] = np.array([10, 20, 20, 50, 100, 200, 300, 300, 500, 500, 500, 500])




    def getOs(self):
        return self.os


    def setBounds(self, imax, jmax ,kmax, fullDepth, t, os):
        bounds = {}
        # eLeft = 0.5*np.ones(jmax)
        # eRight = -0.5*np.ones(jmax)
        # eBottom = np.ones(imax)
        # eTop = np.ones(imax)
        # for i in range(0, imax):
        #     eBottom[i] = 0.5 - i/imax
        #     eTop[i] = 0.5 - i / imax
        # bounds['E'] = [eLeft, eBottom, eRight, eTop ]


        #bounds['U'] = [np.zeros((fullDims[1], os.kmax)), np.zeros((fullDims[0]-1, os.kmax)),\
        #               np.zeros((fullDims[1], os.kmax)), np.zeros((fullDims[0]-1, os.kmax))]
        #bounds['V'] = [np.zeros((fullDims[1]-1, os.kmax)), np.zeros((fullDims[0], os.kmax)), \
        #               np.zeros((fullDims[1]-1, os.kmax)), np.zeros((fullDims[0], os.kmax))]
        return bounds