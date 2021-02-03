import mpi

# # Set wind speeds for our slice based on the current scenario.
# def setAtmo(scenario, fullDims, slice, t, os):
#     windU, windV = scenario.setAtmo(fullDims[0], fullDims[1], t, os)
#     if not windU is None:
#         os.windU[...] = windU[slice[0]:slice[1]-1,slice[2]:slice[3]]
#         os.windV[...] = windV[slice[0]:slice[1], slice[2]:slice[3]-1]


# Set wind speeds for our slice based on the current scenario.
def setAtmo(scenario, comm, fullDims, doMpi, rank, pos, splits, slice, t, os):
    if not doMpi:
        windU, windV = scenario.setAtmo(fullDims[0], fullDims[1], t, os)
        if not windU is None:
            os.windU[...] = windU[...]#[slice[0]:slice[1]-1,slice[2]:slice[3]]
            os.windV[...] = windV[...]#[slice[0]:slice[1], slice[2]:slice[3]-1]
    else:
        if rank==0:
            windUAll, windVAll = scenario.setAtmo(fullDims[0], fullDims[1], t, os)
            # Set status according to whether we have a new atmo field to set:
            if windUAll is None:
                status = 0
            else:
                status = 1
        else:
            status = None

        # Broadcast the status from rank 0 to the others:
        status = mpi.broadcastStatus(comm, rank, status)

        # Scatter the atmo data from rank 0 if we have data:
        if status>0:
            if rank==0:
                windU = mpi.scatter2DField(comm, fullDims, windUAll, rank, pos, splits, slice, [1, 0])
                windV = mpi.scatter2DField(comm, fullDims, windVAll, rank, pos, splits, slice, [0, 1])
            else:
                windU = mpi.scatter2DField(comm, fullDims, None, rank, pos, splits, slice, [1, 0])
                windV = mpi.scatter2DField(comm, fullDims, None, rank, pos, splits, slice, [0, 1])

            os.windU[...] = windU[...]
            os.windV[...] = windV[...]
