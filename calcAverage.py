import oceanState
import mpi
import netcdfStorage

class Average:

    # Initialize object. The OceanState given as argument is used only to
    # get the dimensions of the model domain.
    def __init__(self, os):
        self.avg = oceanState.OceanState(os.imax, os.jmax, os.kmax, includeAncillaries=False)
        self.avg.dz = os.dz.copy()
        self.avg.layerDepths = os.layerDepths.copy()
        self.avg.depth[...] = os.depth[...]
        self.nValues = 0 # The number of times we have recorded values

    def record(self, os):
        self.nValues = self.nValues + 1
        self.avg.U[...] = self.avg.U[...] + os.U[...]
        self.avg.V[...] = self.avg.V[...] + os.V[...]
        self.avg.W[...] = self.avg.W[...] + os.W[...]
        self.avg.T[...] = self.avg.T[...] + os.T[...]
        self.avg.S[...] = self.avg.S[...] + os.S[...]
        self.avg.E[...] = self.avg.E[...] + os.E[...]
        self.avg.X[...] = self.avg.X[...] + os.X[...]

    def doAverage(self):
        self.avg.U = self.avg.U / self.nValues
        self.avg.V = self.avg.V / self.nValues
        self.avg.W = self.avg.W / self.nValues
        self.avg.T = self.avg.T / self.nValues
        self.avg.S = self.avg.S / self.nValues
        self.avg.E = self.avg.E / self.nValues
        self.avg.X = self.avg.X / self.nValues

    def collectAndSaveAverage(self, firstSave, t, doMpi, comm, rank, pos,
                              splits, myHalo, fullDims, fullDepth, sp):
        if rank==0:
            print("Saving average values")

        # All ranks should calculate the average values before we collect all:
        self.doAverage()

        # Before saving, collect full model average. All ranks need to participate here:
        if doMpi:
            osAll = mpi.collectAll(comm, rank, pos, splits, myHalo, self.avg, fullDims, fullDepth, sp)
        else:
            osAll = self.avg

        # The actual save is performed by rank 0 only:
        if rank == 0:

            # Initialize file if this is first save:
            if firstSave:
                firstSave = False
                netcdfStorage.initSaveFile(sp.saveAvgFile, osAll.imax, osAll.jmax, osAll.kmax, osAll.depth,
                                           osAll.layerDepths)
            # Save state:
            netcdfStorage.saveState(sp.saveAvgFile, t, osAll)

        # After saving, we clear out all values in preparation for next round:
        self.nValues = 0
        self.avg.U = 0 * self.avg.U
        self.avg.V = 0 * self.avg.V
        self.avg.W = 0 * self.avg.W
        self.avg.T = 0 * self.avg.T
        self.avg.S = 0 * self.avg.S
        self.avg.E = 0 * self.avg.E
        self.avg.X = 0 * self.avg.X