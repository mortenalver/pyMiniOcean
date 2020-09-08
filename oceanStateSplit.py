from oceanState import OceanState
import numpy as np
import math
import matplotlib.pyplot as plt

# This subclass of OceanState adds extra members that are needed for a
# mode-splitting integration scheme.
class OceanStateSplit(OceanState):

    def __init__(self, imax, jmax, kmax, includeAncillaries=True):
        super(OceanStateSplit, self).__init__(imax, jmax, kmax, includeAncillaries)

        # Define extra cell height matrices at U and V grid locations.
        # DW/DS contain cell heights without elevation, and DWD/DSD
        # will be updated to include elevation.
        self.DW = np.zeros((imax-1,jmax,kmax))
        self.DS = np.zeros((imax,jmax-1,kmax))
        self.DWD = np.zeros((imax - 1, jmax, kmax))
        self.DSD = np.zeros((imax, jmax - 1, kmax))
        self.DWA = np.zeros((imax - 1, jmax))
        self.DSA = np.zeros((imax, jmax - 1))


        self.huSqr = np.zeros((imax-1, jmax))
        self.hvSqr = np.zeros((imax, jmax-1))

        # 2D states for depth times average U/V:
        self.HUA = np.zeros((imax - 1, jmax))
        self.HVA = np.zeros((imax, jmax-1))
        # 2D states for average U/V:
        self.UA = np.zeros((imax - 1, jmax))
        self.VA = np.zeros((imax, jmax-1))
        # 3D states for U/V deviations from average:
        self.UB = np.zeros((imax-1, jmax, kmax))
        self.VB = np.zeros((imax, jmax-1, kmax))

    def calcKmmDzz(self):
        super().calcKmmDzz()

        self.DW[:,:,:] = np.minimum(self.dzz[1:,:,:], self.dzz[0:-1,:,:])
        self.DS[:,:,:] = np.minimum(self.dzz[:,1:,:], self.dzz[:,0:-1,:])
        self.DWD[...] = self.DW[...]
        self.DSD[...] = self.DS[...]

        # Set huSqr:
        for i in range(0, self.imax-1):
            for j in range(0, self.jmax):
                if self.maskU[i,j,0]:
                    self.huSqr[i,j] = math.sqrt(min(self.depth[i,j], self.depth[i+1,j])*9.81)
                self.DWA[i,j] = np.sum(self.DW[i,j,:])
         # Set hvSqr:
        for i in range(0, self.imax):
            for j in range(0, self.jmax-1):
                if self.maskV[i,j,0]:
                    self.hvSqr[i,j] = math.sqrt(min(self.depth[i,j], self.depth[i,j+1])*9.81)
                self.DSA[i, j] = np.sum(self.DS[i, j, :])



    # Based on U and V, calculate the depth-averaged values and the deviations:
    def calc2D3D(self):
        self.HUA[...] = 0
        self.HVA[...] = 0
        for i in range(0, self.imax-1):
            for j in range(0, self.jmax):
                sumD=0
                k=0
                while k<self.kmax and self.maskU[i,j,k]:
                    self.HUA[i,j] = self.HUA[i,j] + self.DWD[i,j,k]*self.U[i,j,k]
                    sumD = sumD + self.DWD[i,j,k]
                    k = k+1
                if k>0:
                    self.UA[i,j] = self.HUA[i,j]/sumD


        for i in range(0, self.imax-1):
            for j in range(0, self.jmax):
                k=0
                while k<self.kmax and self.maskU[i,j,k]:
                    self.UB[i,j,k] = self.U[i,j,k] - self.UA[i,j]
                    k=k+1

        for i in range(0, self.imax):
            for j in range(0, self.jmax-1):
                sumD = 0
                k=0
                while k<self.kmax and self.maskV[i,j,k]:
                    self.HVA[i,j] = self.HVA[i,j] + self.DSD[i,j,k]*self.V[i,j,k]
                    sumD = sumD + self.DSD[i,j,k]
                    k = k+1
                if k>0:
                    self.VA[i,j] = self.HVA[i,j]/sumD

        for i in range(0, self.imax):
            for j in range(0, self.jmax-1):
                k=0
                while k<self.kmax and self.maskV[i,j,k]:
                    self.VB[i,j,k] = self.V[i,j,k] - self.VA[i,j]
                    k=k+1

    def updateCellHeights(self):
        super().updateCellHeights()

        for i in range(0,self.imax-1):
            for j in range(0,self.jmax):
                if self.DW[i,j,0] > 0:
                    self.DWD[i,j,0] = self.DW[i,j,0] + 0.5*(self.E[i,j]+self.E[i+1,j])

        for i in range(0,self.imax):
            for j in range(0,self.jmax-1):
                if self.DS[i,j,0] > 0:
                    self.DSD[i,j,0] = self.DS[i,j,0] + 0.5*(self.E[i,j]+self.E[i,j+1])




    def reduceToSlice(self, slice):
        super().reduceToSlice(slice)

        self.DW = self.DW[slice[0]:slice[1] - 1, slice[2]:slice[3], ...].copy()
        self.DS = self.DS[slice[0]:slice[1], slice[2]:slice[3]-1, ...].copy()
        self.DWD = self.DWD[slice[0]:slice[1] - 1, slice[2]:slice[3], ...].copy()
        self.DSD = self.DSD[slice[0]:slice[1], slice[2]:slice[3]-1, ...].copy()
        self.DWA = self.DWA[slice[0]:slice[1] - 1, slice[2]:slice[3], ...].copy()
        self.DSA = self.DSA[slice[0]:slice[1], slice[2]:slice[3]-1, ...].copy()

        self.huSqr = self.huSqr[slice[0]:slice[1]-1, slice[2]:slice[3], ...].copy()
        self.hvSqr = self.hvSqr[slice[0]:slice[1], slice[2]:slice[3]-1, ...].copy()

        self.HUA = self.HUA[slice[0]:slice[1] - 1, slice[2]:slice[3], ...].copy()
        self.HVA = self.HVA[slice[0]:slice[1], slice[2]:slice[3]-1, ...].copy()
        self.UA = self.UA[slice[0]:slice[1] - 1, slice[2]:slice[3], ...].copy()
        self.VA = self.VA[slice[0]:slice[1], slice[2]:slice[3]-1, ...].copy()
        self.UB = self.UB[slice[0]:slice[1] - 1, slice[2]:slice[3], ...].copy()
        self.VB = self.VB[slice[0]:slice[1], slice[2]:slice[3] - 1, ...].copy()

