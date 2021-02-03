import netcdfStorage
from oceanState import *
from oceanStateSplit import *
from simSettings import *
from verticalSpeeds import *
from setBounds import *
import openArea, smallScale, upwelling, real, fjord
import atmo, freshwater, advection, integration, calcAverage, split
import time
from mpi4py import MPI
import mpi
import utils
import sys

#scenario = openArea.OpenArea()
#scenario = smallScale.SmallScale()
#scenario = upwelling.Upwelling()
scenario = real.Real()
#scenario = fjord.Fjord()
coldStart = True

plotInt = -1 # Interval (samples) between updating plot(s). Set to -1 to disable plotting.

comm = MPI.COMM_WORLD

pos = (0,0)
splits = (1,1)
# Discover the MPI settings, our own rank, and our position in the tile grid:
#mpisize, rank, splits, pos = mpi.mpiSettings(comm)
#doMpi = mpisize > 1
#if doMpi:
#    print("MPI: rank="+str(rank)+", pos="+str(pos)+", splits="+str(splits))
#    #if rank==0:
#    #    statusfile = open('status.txt', 'w')
#else:
#    print("One process only, disabling MPI communication.")

# Initialize sim settings:
sp = SimSettings()

if coldStart:
    scenario.initialize(sp) # Set up scenario
    os = scenario.getOs() # Get OceanState object from the scenario
    utils.adaptDepthField(os) # Make some adjustments to ensure the depth field works well with bounds
    if sp.passiveTracer:
        scenario.initPassiveTracer()


else:
    # Load saved model state:
    scenario.initialize(sp) # Set up scenario
    os = netcdfStorage.loadState(sp, 'f:/work/pyMiniOcean/fjord1.nc', 40) # Load OceanState from file
    utils.adaptDepthField(os)  # Make some adjustments to ensure the depth field works well with bounds
    scenario.os = os # Set the scenario's OceanState object
    #os.U = 0*os.U
    #os.V = 0*os.V


# If we are using mode splitting, after initialization we need to calculate the initial
# values of the depth integrated speeds and deviations:
if sp.modeSplittingOn:
    os.calc2D3D()

    #plt.figure()
    #plt.pcolor(np.transpose(os.UB[:,:,0])), plt.colorbar()
    #plt.show()

computeVerticalSpeeds(os, sp) # Calculate vertical speeds based on initial horizontal speeds

t0 = time.time() # For timing purposes only


# if doMpi:
#     # Make note of the full model dimensions, then determine which slice this instance
#     # is going to cover.
#     fullDims = (os.imax, os.jmax, os.kmax) # Dimensions of full model domain
#     fullDepth = os.depth
#     sliceNoHalo = mpi.getSlice(splits, pos, os.imax, os.jmax)
#     slice = mpi.getSliceWithHalo(splits, pos, os.imax, os.jmax)
#     myHalo = (sliceNoHalo[0]-slice[0], slice[1]-sliceNoHalo[1], \
#               sliceNoHalo[2] - slice[2], slice[3] - sliceNoHalo[3])
#     # Slice down the model domain to the correct slice:
#     os.reduceToSlice(slice)
# else:
fullDims = (os.imax, os.jmax, os.kmax)  # Dimensions of full model domain
fullDepth = os.depth
slice = (0, os.imax, 0, os.jmax)
myHalo = None

rank = comm.Get_rank()
myInputFileName = "input_"+(str(rank)).zfill(3)+".nc"
print(myInputFileName)

nSamples = int(sp.tEnd/sp.dt)
saveIntSamples = int(sp.saveIntS/sp.dt)
if rank==0:
    print("Save interval = "+str(saveIntSamples))
saveCount = 0
firstSave = True
calcMixingCount = 0
pltCount = 0
if sp.recordAverages:
    saveAverageCount = 0
    saveAvgIntSamples = int(sp.saveAvgIntS / sp.dt)
    firstAvgSave = True
    aver = calcAverage.Average(os)

#if plotInt > 0:
    # myfig = plt.figure()
    # g1 = myfig.add_subplot(221)
    # p1 = g1.pcolor(os.U[:,:-1,0])
    # g2 = myfig.add_subplot(222)
    # g3 = myfig.add_subplot(223)
    # g4 = myfig.add_subplot(224)

#aU, aV = utils.vertAverageUBVB(os, 20, 20)
#print("aU="+str(aU)+", aV="+str(aV))


for sample in range(0,nSamples):

    t = (sample+1)*sp.dt # Current model time

    # Update cell heights to account for updated elevation:
    os.updateCellHeights()
    # Update densities to account for updated T and S:
    os.updateRho()
    # Update mixing coefficients (every nth sample):
    calcMixingCount = calcMixingCount + 1
    if sample==0 or calcMixingCount == sp.calcMixingInterval:
        calcMixingCount = 0
        os.calcVerticalMixingRichardson(sp)
        os.calcHorizontalDiffusivitySmagorinsky(sp)
        if sp.modeSplittingOn:
            os.calcHorizontalDiffusivitySmagorinsky2D(sp)


    # Set boundary values:
    setBounds(scenario, fullDims, fullDepth, pos, splits, slice, t, os, sp)

    # Set freshwater:
    if sp.freshwaterOn:
        freshwater.addFreshwater(scenario, fullDims, pos, splits, slice, t, os, sp)

    # Set atmo values:
    if sp.atmoOn:
        atmo.setAtmo(scenario, fullDims, pos, splits, slice, t, os)

    # Calculate T and S for next time step:
    print(os.T[1, 17:20, 0])
    print(os.T_next[1, 17:20, 0])
    advection.advectTempSalt(os, sp)
    print(os.V[1, 17:20, 0])
    print(os.T[1, 17:20, 0])

    # Advect passive tracer:
    if sp.passiveTracer:
        advection.advectPassiveTracer(os, sp)
        os.X[1:-1,1:-1,:] = os.X_next[1:-1,1:-1,:]

    # Calculate U, V and E for next time step. Calculation method depends on simulation settings:
    if sp.modeSplittingOn:
        # Mode splitting means a split barotropic-baroclinic integration scheme that is more efficient:
        split.integrate(os, sp, scenario, fullDims, fullDepth, pos, splits, slice, t)

        #aU, aV = utils.vertAverageUBVB(os, 20, 20)
        #print("aU=" + str(aU) + ", aV=" + str(aV))
    else:
        # Call the chosen non-split integration scheme:
        if not sp.useRk:
            integration.forwardEuler(os, sp, t)
        else:
            integration.rk4Integration(os, sp, t)

    os.T[1:-1,1:-1,:] = os.T_next[1:-1,1:-1,:]
    os.S[1:-1, 1:-1, :] = os.S_next[1:-1, 1:-1, :]

    ####################### TIME STEP DONE ########################

    # If active, record and possibly save averages:
    if sp.recordAverages:
        saveAverageCount = saveAverageCount+1
        aver.record(os)
        if saveAverageCount == saveAvgIntSamples:
            saveAverageCount = 0
            aver.collectAndSaveAverage(firstAvgSave, t, doMpi, comm, rank, pos, splits,
                                       myHalo, fullDims, fullDepth, sp)
            firstAvgSave = False



    # See if we should save state:
    saveCount = saveCount+1
    if saveCount == saveIntSamples:
        saveCount = 0
        t1 = time.time()
        if rank == 0:
            print("Time = "+str(t/3600)+" h. S since last save = " + str(t1 - t0))
        t0 = t1




