
# Set wind speeds for our slice based on the current scenario.
def setAtmo(scenario, fullDims, pos, splits, slice, t, os):

    windU, windV = scenario.setAtmo(fullDims[0], fullDims[1], t, os)
    if not windU is None:
        os.windU = windU[slice[0]:slice[1]-1,slice[2]:slice[3]]
        os.windV = windV[slice[0]:slice[1], slice[2]:slice[3]-1]
