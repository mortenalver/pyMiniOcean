from mpi4py import MPI
import math
import numpy as np
import oceanState

halo = 2 # The extra grid cells included around each tile representing the
    # interface with neightbouring tiles
dummyOs = oceanState.OceanState(2, 2, 2) # Dummy OceanState object

def mpiSettings(comm):
    rank = comm.Get_rank()
    size = comm.Get_size()

    if size==1:
        splits = (1, 1) # TEST
    elif size==2:
        splits = (2, 1)
    elif size == 3:
        splits = (3, 1)
    elif size==4:
        splits = (2, 2)
    elif size==6:
        splits = (3, 2)
    elif size==8:
        splits = (4, 2)
    elif size==9:
        splits = (3, 3)
    elif size==12:
        splits = (4, 3)
    elif size==16:
        splits = (4, 4)
    elif size == 20:
        splits = (5, 4)
    elif size == 32:
        splits = (8, 4)
    elif size == 40:
        splits = (8, 5)
    x = rank - splits[0]*int(rank/splits[0])
    y = int(rank/splits[0])

    return size, rank, splits, (x, y)

# Compute which slice of the complete area to cover for a given pos(x,y)
# of the tile set dimensions splits (x,y), and total domain size imax,jmax.
def getSlice(splits, pos, imax, jmax):
    #print("pos = "+str(pos))
    istart = int(math.floor(pos[0]*imax/splits[0]))
    if pos[0]>0:
        istart = istart+1
    iend = min(imax, 1+(math.floor((pos[0]+1) * imax / splits[0])))
    #print("istart = " + str(istart))
    #print("iend = " + str(iend))
    jstart = int(math.floor(pos[1] * jmax / splits[1]))
    if pos[1]>0:
        jstart = jstart+1
    jend = min(jmax, 1+int(math.floor((pos[1]+1)*jmax / splits[1])))
    return (istart, iend, jstart, jend)


# Compute which slice of the complete area to cover for a given pos(x,y)
# of the tile set dimensions splits (x,y), and total domain size imax,jmax.
# The slice is given a halo of extra cells along internal borders, so the
# slices overlap. This is to account for the need in each subdomain to know
# about the state right outside its border, and the halos are the areas where
# we need to interchange information using MPI functions.
def getSliceWithHalo(splits, pos, imax, jmax):
    slice = getSlice(splits, pos, imax, jmax)
    return (max(0, slice[0]-halo), min(imax, slice[1]+halo), \
            max(0, slice[2] - halo), min(jmax, slice[3] + halo))


def communicateVar(comm, rank, pos, splits, var, edgOffs):
    # Make note of the shape and number of dimensions of the variable:
    shape = var.shape
    dims = len(shape)
    # Go through the grid of tiles, sending the right edge of subdomains to
    # overwrite the left edge halo of the next door domains:
    for i in range(0,splits[0]-1):
        for j in range(0, splits[1]):
            # (i,j) denotes the position we want to send from.
            # If we are the subdomain that should send data, do it:
            if pos == (i, j):
                data = var[shape[0]+edgOffs[0]-2*halo:shape[0]+edgOffs[0]-halo,...]
                data2 = np.empty(data.shape)
                data2[...] = data[...]
                comm.Send([data2, MPI.FLOAT], dest=rank+1, tag=10)
            # If we are the next subdomain that should receive the data, do it:
            elif pos == (i+1, j):
                if dims==2:
                    data = np.empty((halo, shape[1]))
                else:
                    data = np.empty((halo, shape[1], shape[2]))
                comm.Recv([data, MPI.FLOAT], source=rank-1, tag=10)
                var[0:halo,...] = data
                #print(data)

    # Do the same, now sending the left edge of subdomains to overwrite the
    # right edge halo of the next domains:
    for i in range(1, splits[0]):
        for j in range(0, splits[1]):
            if pos == (i, j):
                data = var[halo-edgOffs[0]:2*halo-edgOffs[0],...]
                data2 = np.empty(data.shape)
                data2[...] = data[...]
                comm.Send([data2, MPI.FLOAT], dest=rank-1, tag=10)
            elif pos == (i-1, j):
                if dims==2:
                    data = np.empty((halo, shape[1]))
                else:
                    data = np.empty((halo, shape[1], shape[2]))
                comm.Recv([data, MPI.FLOAT], source=rank+1, tag=10)
                var[shape[0]-halo:shape[0],...] = data

    # Repeat vertically, sending upwards:
    for i in range(0,splits[0]):
        for j in range(0, splits[1]-1):
            # (i,j) denotes the position we want to send from.
            # If we are the subdomain that should send data, do it:
            if pos == (i, j):
                data = var[:, shape[1]+edgOffs[1]-2*halo:shape[1]+edgOffs[1]-halo, ...]
                data2 = np.empty(data.shape)
                data2[...] = data[...]
                comm.Send([data2, MPI.FLOAT], dest=rank + splits[0], tag=10)
            # If we are the next subdomain that should receive the data, do it:
            elif pos == (i, j+1):
                if dims==2:
                    data = np.empty((shape[0], halo))
                else:
                    data = np.empty((shape[0], halo, shape[2]))
                comm.Recv([data, MPI.FLOAT], source=rank-splits[0], tag=10)
                var[:,0:halo,...] = data
    # Sending downwards:
    for i in range(0,splits[0]):
        for j in range(1, splits[1]):
            # (i,j) denotes the position we want to send from.
            # If we are the subdomain that should send data, do it:
            if pos == (i, j):
                data = var[:, halo-edgOffs[1]:2*halo-edgOffs[1], ...]
                data2 = np.empty(data.shape)
                data2[...] = data[...]
                comm.Send([data2, MPI.FLOAT], dest=rank-splits[0], tag=10)
            # If we are the next subdomain that should receive the data, do it:
            elif pos == (i, j-1):
                if dims==2:
                    data = np.empty((shape[0], halo))
                else:
                    data = np.empty((shape[0], halo, shape[2]))
                comm.Recv([data, MPI.FLOAT], source=rank+splits[0], tag=10)
                var[:,shape[1]-halo:shape[1],...] = data


# Do exchange of information between tiles in order to update
# each tile's halo with correct values from their neighbours:
def communicate(comm, rank, pos, splits, os, sp):
    communicateVar(comm, rank, pos, splits, os.E, (0,0))
    communicateVar(comm, rank, pos, splits, os.U, (1,0))
    communicateVar(comm, rank, pos, splits, os.V, (0, 1))
    communicateVar(comm, rank, pos, splits, os.T, (0,0))
    communicateVar(comm, rank, pos, splits, os.S, (0, 0))
    if sp.modeSplittingOn:
        communicateVar(comm, rank, pos, splits, os.HUA, (1, 0))
        communicateVar(comm, rank, pos, splits, os.HVA, (0, 1))
        communicateVar(comm, rank, pos, splits, os.UA, (1, 0))
        communicateVar(comm, rank, pos, splits, os.VA, (0, 1))
        communicateVar(comm, rank, pos, splits, os.UB, (1, 0))
        communicateVar(comm, rank, pos, splits, os.VB, (0, 1))

# Collector function handling a single variable. If rank==0, this function
# will be called with allVar as the OceanState into which all tiles should
# be collected. If rank>0, allVar will be a dummy OceanState, as this instance
# will only send information to the rank==0 instance.
def collectOneVar(comm, rank, pos, splits, myHalo, fullDims, var, allVar):
    count = -1
    threeDim = len(var.shape)>2
    if threeDim:
        kmax = var.shape[2]
    else:
        kmax = 0

    for j in range(0, splits[1]):
        for i in range(0, splits[0]):
            count = count+1
            slice = getSlice(splits, (i, j), fullDims[0], fullDims[1])
            if rank>0:
                # If this is my position, I should send:
                if i==pos[0] and j==pos[1]:
                    #print("rank=" + str(rank) + ", slice=" + str(slice))
                    #print("pos="+str((i,j))+", rank "+str(rank)+" sends. myHalo="+str(myHalo))
                    data = var[myHalo[0]:myHalo[0]+slice[1]-slice[0], myHalo[2]:myHalo[2]+slice[3]-slice[2], ...]
                    #print("rank=" + str(rank)+", count="+str(count)+  ", griddim=" + str(var.shape)+", data="+str(data.shape))
                    data2 = np.empty(data.shape)
                    data2[...] = data[...]
                    comm.Send([data2, MPI.FLOAT], dest=0, tag=count)

            else:
                if count > 0:
                    #print("rank 0, allVar="+str(allVar.shape)+", slice="+str(slice))
                    if threeDim:
                        data = np.empty((min(allVar.shape[0],slice[1]) - slice[0], min(allVar.shape[1],slice[3]) - slice[2], kmax))
                    else:
                        data = np.empty((slice[1]-slice[0], slice[3]-slice[2]))
                    #print("pos=" + str((i, j)) + ", rank " + str(rank) + " receives from "+str(count)+". slice=" + str(slice)+", datadim="+str(data.shape))
                    comm.Recv([data, MPI.FLOAT], source=count, tag=count)
                    print(data.shape)
                    print((allVar[slice[0]:slice[0]+data.shape[0],slice[2]:slice[2]+data.shape[1],...]).shape)
                    #print("rank="+str(rank)+", received data="+str(data.shape))
                    #if i!=1 or j!=1:
                    allVar[slice[0]:slice[0]+data.shape[0],slice[2]:slice[2]+data.shape[1],...] = data
                else:
                    #print("rank 0 should update my own slice")
                    #allVar[slice[0]:slice[1],slice[2]:slice[3],...] = var[myHalo[0]:var.shape[0]-myHalo[1],myHalo[2]:var.shape[1]-myHalo[3],...]
                    allVar[slice[0]:slice[1], slice[2]:slice[3], ...] = var[myHalo[0]:myHalo[0]+slice[1]-slice[0],
                                                                        myHalo[2]:myHalo[2]+slice[3]-slice[2], ...]



# If rank==0, build a full model state based on collecting information from
# all tiles. If rank>0, supply the needed information.
def collectAll(comm, rank, pos, splits, myHalo, os, fullDims, fullDepth, sp):
    if rank == 0:
        osAll = oceanState.OceanState(fullDims[0], fullDims[1], fullDims[2])
        osAll.depth = fullDepth
        osAll.layerDepths = os.layerDepths
    else:
        osAll = dummyOs

    collectOneVar(comm, rank, pos, splits, myHalo, fullDims, os.E, osAll.E)
    collectOneVar(comm, rank, pos, splits, myHalo, fullDims, os.U, osAll.U)
    collectOneVar(comm, rank, pos, splits, myHalo, fullDims, os.V, osAll.V)
    collectOneVar(comm, rank, pos, splits, myHalo, fullDims, os.W, osAll.W)
    collectOneVar(comm, rank, pos, splits, myHalo, fullDims, os.T, osAll.T)
    collectOneVar(comm, rank, pos, splits, myHalo, fullDims, os.S, osAll.S)
    #collectOneVar(comm, rank, pos, splits, myHalo, fullDims, os.windU, osAll.windU)
    #collectOneVar(comm, rank, pos, splits, myHalo, fullDims, os.windV, osAll.windV)

    if sp.passiveTracer:
        collectOneVar(comm, rank, pos, splits, myHalo, fullDims, os.X, osAll.X)

    return osAll