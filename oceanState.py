import numpy as np
import math
from utils import getUV, getUV2
import matplotlib.pyplot as plt

class OceanState:

    minBotLayerDepth = 1.0

    def __init__(self, imax, jmax, kmax, includeAncillaries=True):
        self.includeAncillaries = includeAncillaries
        self.imax = imax
        self.jmax = jmax
        self.kmax = kmax
        self.depth = np.zeros((imax, jmax))
        self.dz = np.zeros((kmax))
        self.layerDepths = np.zeros((kmax))
        self.midLayerDepths = np.zeros((kmax))
        self.U = np.zeros((imax-1,jmax,kmax))
        self.V = np.zeros((imax,jmax-1,kmax))
        self.W = np.zeros((imax,jmax,kmax+1))
        self.T = np.zeros((imax,jmax,kmax))
        self.S = np.zeros((imax,jmax,kmax))
        self.E = np.zeros((imax,jmax))
        self.cellHeights = np.zeros((imax,jmax,kmax))
        self.kmm = np.zeros((imax, jmax), int)
        self.X = np.zeros((imax,jmax,kmax)) # Passive tracer
        if includeAncillaries:
            self.U_next = np.zeros((imax - 1, jmax, kmax))
            self.maskU = np.zeros((imax - 1, jmax, kmax), int)
            self.V_next = np.zeros((imax, jmax - 1, kmax))
            self.maskV = np.zeros((imax, jmax - 1, kmax), int)
            self.T_next = np.zeros((imax, jmax, kmax))
            self.S_next = np.zeros((imax, jmax, kmax))
            self.E_next = np.zeros((imax, jmax))
            self.windU = np.zeros((imax-1, jmax))
            self.windV = np.zeros((imax, jmax-1))
            self.K_v = np.zeros((imax, jmax, kmax-1))
            self.AH = np.zeros((imax, jmax, kmax))
            self.AM2D = np.zeros((imax, jmax))
            self.X_next = np.zeros((imax, jmax, kmax))  # Passive tracer
            self.dzz = np.zeros((imax,jmax,kmax))
            self.rho = np.zeros((imax, jmax, kmax))


    # Slice all distributed states down to a given slice:
    def reduceToSlice(self, slice):
        self.imax = slice[1]-slice[0]
        self.jmax = slice[3]-slice[2]
        self.depth = self.depth[slice[0]:slice[1], slice[2]:slice[3]].copy()
        self.U = self.U[slice[0]:slice[1]-1, slice[2]:slice[3], ...].copy()
        self.V = self.V[slice[0]:slice[1], slice[2]:slice[3]-1, :].copy()
        self.W = self.W[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
        self.T = self.T[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
        self.S = self.S[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
        self.E = self.E[slice[0]:slice[1], slice[2]:slice[3]].copy()
        self.cellHeights = self.cellHeights[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
        self.kmm = self.kmm[slice[0]:slice[1], slice[2]:slice[3]].copy()
        self.X = self.X[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
        if self.includeAncillaries:
            self.U_next = self.U_next[slice[0]:slice[1] - 1, slice[2]:slice[3], ...].copy()
            self.maskU = self.maskU[slice[0]:slice[1]-1, slice[2]:slice[3], :].copy()
            self.V_next = self.V_next[slice[0]:slice[1], slice[2]:slice[3] - 1, :].copy()
            self.maskV = self.maskV[slice[0]:slice[1], slice[2]:slice[3] - 1, :].copy()
            self.T_next = self.T_next[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
            self.S_next = self.S_next[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
            self.E_next = self.E_next[slice[0]:slice[1], slice[2]:slice[3]].copy()
            self.windU = self.windU[slice[0]:slice[1]-1, slice[2]:slice[3]].copy()
            self.windV = self.windV[slice[0]:slice[1], slice[2]:slice[3]-1].copy()
            self.K_v = self.K_v[slice[0]:slice[1], slice[2]:slice[3]].copy()
            self.AH = self.AH[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
            self.AM2D = self.AM2D[slice[0]:slice[1], slice[2]:slice[3]].copy()
            self.X_next = self.X_next[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
            self.dzz = self.dzz[slice[0]:slice[1], slice[2]:slice[3], ...].copy()
            self.rho = self.rho[slice[0]:slice[1], slice[2]:slice[3], ...].copy()

    def calcKmmDzz(self):
        # Calculate layerDepths and midLayerDepths:
        sum = 0
        for k in range(0,self.kmax):
            sum = sum+self.dz[k]
            self.layerDepths[k] = sum
            self.midLayerDepths[k] = sum - 0.5*self.dz[k]


        # Calculate kmm, dzz and cellHeights:
        self.kmm[:,:] = 0
        for i in range(0,self.imax):
            for j in range(0,self.jmax):
                if self.depth[i,j]>0:
                    if self.depth[i,j] < self.layerDepths[0]:
                        self.depth[i,j] = self.dz[0]
                        self.dzz[i,j,0] = self.dz[0]
                        self.kmm[i,j] = 1
                    else:
                        for k in range(0,self.kmax):
                            if self.depth[i,j] > self.layerDepths[k]:
                                self.dzz[i,j,k] = self.dz[k]
                            else:
                                self.dzz[i, j, k] = self.depth[i,j] - self.layerDepths[k-1]
                                if self.dzz[i, j, k] < self.minBotLayerDepth:
                                    self.depth[i,j] = self.layerDepths[k-1] + self.minBotLayerDepth
                                    self.dzz[i,j,k] = self.minBotLayerDepth
                                self.kmm[i,j] = k+1
                                break
                        if self.depth[i,j] > 0 and self.kmm[i,j]==-1:
                            self.depth[i,j] = self.layerDepths[self.kmax-1]
                            self.kmm[i,j] = self.kmax


        self.cellHeights = self.dzz.copy()

        self.maskU[...] = 0
        self.maskV[...] = 0

        for i in range(0,self.imax-1):
            for j in range(0,self.jmax):
                for k in range(0,self.kmax):
                    if self.kmm[i,j] > k and self.kmm[i+1,j] > k:
                        self.maskU[i,j,k] = 1
                    else:
                        self.U[i,j,k] = math.nan
                    #if i<=1 and j==54 and k==5:
                    #    print(str(i)+", kmm="+str(self.kmm[i,j])+", maskU="+str(self.maskU[i,j,k])+", U="+str(self.U[i,j,k]))

        for i in range(0,self.imax):
            for j in range(0,self.jmax-1):
                for k in range(0,self.kmax):
                    if self.kmm[i,j] > k and self.kmm[i,j+1] > k:
                        self.maskV[i,j,k] = 1
                    else:
                        self.V[i,j,k] = math.nan


    def updateCellHeights(self):
        for i in range(0,self.imax):
            for j in range(0,self.jmax):
                if self.kmm[i,j] >= 0:
                    if not math.isnan(self.E[i,j]):
                        self.cellHeights[i,j,0] = self.dzz[i,j,0] + self.E[i,j]



    def dens(self, S, T):
        T3 = T*T*T
        return 999.842594 + 6.793952E-2 * T \
             - 9.095290E-3 * T * T + 1.001685E-4 * T3 \
             - 1.120083E-6 * T3 * T + 6.536332E-9 * T3 * T * T \
             + 8.24493E-1 * S - 4.0899E-3 * T * S \
             + 7.6438E-5 * T * T * S - 8.2467E-7 * T3 * S \
             + 5.3875E-9 * T3 * T * S - 5.72466E-3 * (math.pow(S, 1.5)) \
             + 1.0227E-4 * T * (math.pow(S, 1.5)) - 1.6546E-6 * T * T * (math.pow(S, 1.5)) \
             + 4.8314E-4 * S * S

    def updateRho(self):
        for i in range(0,self.imax):
            for j in range(0,self.jmax):
                for k in range(0,self.kmm[i,j]):
                    self.rho[i,j,k] = self.dens(self.S[i,j,k], self.T[i,j,k])

    # Compute vertical mixing coefficients using the Richardson number scheme
    # (see e.g. Sundfjord et al. 2008; Vertical mixing in the marginal ice
    # zone of the northern Barents Sea—Results from numerical model
    # experiments).
    # Coefficients are defined in the cell centre, horizontally, and on
    # the edges, vertically.
    # TODO: skulle sjekket om Richardson-tallet beregnes riktig. Får alltid 0 eller maks mixing???
    def calcVerticalMixingRichardson(self, sp):
        for i in range(0,self.imax):
            for j in range(0,self.jmax):
                lowerDepth = 0
                for k in range(0,self.kmax-1):
                    lowerDepth = lowerDepth + self.cellHeights[i,j,k]
                    centerDepth = lowerDepth - 0.5*self.cellHeights[i,j,k]
                    if k<self.kmm[i,j]:
                        k_w = 0.028*((sp.H_wave*sp.H_wave)/sp.T_wave)*math.exp(-0.8*centerDepth/sp.H_wave)
                        meanCellHeights = 0.5*(self.cellHeights[i,j,k] + self.cellHeights[i,j,k+1])
                        # Calculate the vertical gradient of the
                        # density. If density is increasing downwards
                        # we want a positive value:
                        d_rho_dz = (self.rho[i,j,k+1] - self.rho[i,j,k])/meanCellHeights
                        meanU_above = 0.5*(getUV(self.U,i-1,j,k,0) + getUV(self.U,i,j,k,0))
                        meanU_below = 0.5*(getUV(self.U,i-1,j,k+1,0) + getUV(self.U,i,j,k+1,0))
                        meanV_above = 0.5*(getUV(self.V,i,j-1,k,0) + getUV(self.V,i,j,k,0))
                        meanV_below = 0.5*(getUV(self.V,i,j-1,k+1,0) + getUV(self.V,i,j,k+1,0))
                        d_U_dz2 = (math.pow(meanU_above - meanU_below, 2) + math.pow(meanV_above - meanV_below,2))/ \
                                   (meanCellHeights*meanCellHeights)
                        d_U_dz2 = max(d_U_dz2, 1e-6)
                        #if d_U_dz2 == 0:
                        #    Ri = 0
                        #else:
                        Ri = (9.81/self.rho[i,j,k])*d_rho_dz/d_U_dz2
                        if math.isnan(Ri):
                            Ri = 0

                        self.K_v[i,j,k] = sp.KVm*(math.atan(sp.G_vmix*(sp.Ri0-Ri))/math.pi + 0.5) + k_w


                    else:
                        self.K_v[i,j,k] = 0

    # Calculate horizontal mixing coefficients according to Smagorinsky (...)
    def calcHorizontalDiffusivitySmagorinsky(self, sp):
        dx2 = sp.dx*sp.dx
        amLim= 0.0625*dx2/(2*sp.dt)
        for i in range(1,self.imax-1):
            for j in range(1, self.jmax-1):
                for k in range(0, self.kmm[i,j]):
                    AM = sp.CM*dx2*math.sqrt(
                        math.pow(getUV(self.U,i,j,k,0)-getUV(self.U,i-1,j,k,0)/sp.dx,2)
                        + math.pow((getUV(self.V,i,j,k,0)-getUV(self.V,i,j-1,k,0))/sp.dx,2)
                        +.5*math.pow(.25*(getUV(self.U,i-1,j-1,k,0) + getUV(self.U,i-1,j+1,k,0)
                        - getUV(self.U,i,j-1,k,0) - getUV(self.U,i,j+1,k,0))/sp.dx
                        + .25*(getUV(self.V,i-1,j-1,k,0) + getUV(self.V,i+1,j-1,k,0)- getUV(self.V,i-1,j,k,0)
                            - getUV(self.V,i+1,j,k,0))/sp.dx, 2))

                    self.AH[i,j,k] = min((sp.CH/sp.CM)*AM, amLim)

        self.AH[0,:,:] = self.AH[1,:,:]
        self.AH[-1,:,:] = self.AH[2,:,:]
        self.AH[:,0,:] = self.AH[:,1,:]
        self.AH[:,-1,:] = self.AH[:,-2,:]
        self.AH[0,0,:] = self.AH[1,1,:]
        self.AH[0,-1,:] = self.AH[1,-2,:]
        self.AH[-1,0,:] = self.AH[-2,1,:]
        self.AH[-1,-1,:] = self.AH[-2,-2,:]

    # Calculate 2D horizontal mixing coefficients according to Smagorinsky (...)
    def calcHorizontalDiffusivitySmagorinsky2D(self, sp):
        dx2 = sp.dx * sp.dx
        amLim= 0.0625*dx2/(2*sp.dt)
        for i in range(1,self.imax-1):
            for j in range(1, self.jmax-1):
                AM = sp.CM2D*dx2*math.sqrt(
                    math.pow(getUV2(self.UA,i,j,0)-getUV2(self.UA,i-1,j,0)/sp.dx,2)
                    + math.pow((getUV2(self.VA,i,j,0)-getUV2(self.VA,i,j-1,0))/sp.dx,2)
                    +.5*math.pow(.25*(getUV2(self.UA,i-1,j-1,0) + getUV2(self.UA,i-1,j+1,0)
                    - getUV2(self.UA,i,j-1,0) - getUV2(self.UA,i,j+1,0))/sp.dx
                    + .25*(getUV2(self.VA,i-1,j-1,0) + getUV2(self.VA,i+1,j-1,0)- getUV2(self.VA,i-1,j,0)
                        - getUV2(self.VA,i+1,j,0))/sp.dx, 2))

                self.AM2D[i,j] = min(AM, amLim)

        self.AM2D[0,:] = self.AM2D[1,:]
        self.AM2D[-1,:] = self.AM2D[2,:]
        self.AM2D[:,0] = self.AM2D[:,1]
        self.AM2D[:,-1] = self.AM2D[:,-2]
        self.AM2D[0,0] = self.AM2D[1,1]
        self.AM2D[0,-1] = self.AM2D[1,-2]
        self.AM2D[-1,0] = self.AM2D[-2,1]
        self.AM2D[-1,-1] = self.AM2D[-2,-2]

