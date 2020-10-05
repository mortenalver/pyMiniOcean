from scenario import *
from oceanState import *
import utils
import math
import netcdfStorage
import datetime

class Real(Scenario):

    os = None

    sinmodTile = [75,105, 35, 60]
    sinmodInitFile = 'C:/SINMOD_operative/gin/gin_PhysStates.nc'
    sinmodAtmoFile = 'C:/SINMOD_operative/gin/gin_atmo_converted.nc'
    initTime = None
    atmoTimes = None
    atmoTimeStep = None
    atmoSample = -1

    def initialize(self, sp):
        sp.dt = 1200
        sp.nsub = 25
        sp.dx = 20000
        sp.KBi = 1.5e12  # Biharmonic constant

        self.E_floating = False
        self.U_floating = True
        self.V_floating = True

        #self.os = netcdfStorage.loadState(sp, 'data/real_init.nc', 0)
        #self.os = netcdfStorage.loadSINMODState(sp, 'C:/SINMOD_operative/gin/gin_PhysStates.nc', 0, [25, 85, 25, 80])
        self.os, self.initTime = netcdfStorage.loadSINMODState(sp, self.sinmodInitFile, 75, self.sinmodTile)
        print("Initializing at time: "+str(self.initTime))
        #self.os = OceanState(300, 235, 5)
        #self.os.dz[:] = np.array([10, 20, 20, 50, 100, 200, 300, 300, 500, 500, 500, 500])

        self.atmoTimes, self.atmoTimeStep = netcdfStorage.getSINMODAtmoTimes(self.sinmodAtmoFile)
        print(self.atmoTimes)



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

    # The setAtmo method gets dimensions for the full domain as input, and should return
    # matrices containing the U and V components of the wind field.
    def setAtmo(self, imax, jmax, t, os):
        roundedT = round(self.atmoTimeStep*round(t/self.atmoTimeStep))
        modelTime = self.initTime + datetime.timedelta(seconds=roundedT)
        tIndex = self.atmoTimes.index(modelTime)
        if self.atmoSample==-1 or tIndex>self.atmoSample:
            self.atmoSample = tIndex
            WU, WV = netcdfStorage.loadSINMODAtmo(self.sinmodAtmoFile, self.atmoSample, self.sinmodTile)
            return WU, WV
        else:
            return None, None
