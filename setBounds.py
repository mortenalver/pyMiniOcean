import numpy as np

def copyInnerValues(pos, splits, var):
    # Left:
    if pos[0] == 0:
        var[0, :, ...] = var[1, :, ...]
    # Bottom:
    if pos[1] == 0:
        var[:, 0, ...] = var[:, 1, ...]
    # Right:
    if pos[0] == splits[0] - 1:
        var[var.shape[0] - 1, :, ...] = var[var.shape[0] - 2, :, ...]
    # Top:
    if pos[1] == splits[1] - 1:
        var[:, var.shape[1] - 1, ...] = var[:, var.shape[1] - 2, ...]


# Apply the bounds for a given variable retrieved from the scenario (for the full domain)
# to our slice of the domain.
def setBoundsVar(bounds, fullDims, pos, splits, slice, varName, var, floating, offset):
    # First check if bounds have been specified for this variable:
    if not varName in bounds:
        if not floating:
            return
        else:
            # No bounds given, but "floating" is set to true. This means we should
            # copy the values right inside the edges out to the edges:
            copyInnerValues(pos, splits, var)
            return

    bnd = bounds[varName] # [left, bottom, right, upper]

    # Left:
    if pos[0]==0:
        var[0,:,...] = bnd[0][slice[2]:slice[3]-offset[1],...]
    # Bottom:
    if pos[1]==0:
        var[:,0,...] = bnd[1][slice[0]:slice[1]-offset[0],...]
    # Right:
    if pos[0]==splits[0]-1:
        var[var.shape[0]-1,:,...] = bnd[2][slice[2]:slice[3]-offset[1],...]
    # Top:
    if pos[1]==splits[1]-1:
        var[:,var.shape[1]-1,...] = bnd[3][slice[0]:slice[1]-offset[0],...]




# Update boundary values, depending on the current scenario. Because of the MPI parallelization, we
# first need to let the scenario set boundaries as if we were working on the full model domain, then
# we pick out the parts that we need.
# If we are not adjacent to any boundaries, we skip the process.
def setBounds(scenario, fullDims, fullDepth, pos, splits, slice, t, os, sp):
    if pos[0]>0 and pos[0]<splits[0]-1 and pos[1]>0 and pos[1]<splits[1]-1:
        return # This slice does not have any external boundaries
    bounds = scenario.setBounds(fullDims[0], fullDims[1], fullDims[2], fullDepth, t, os)

    setBoundsVar(bounds, fullDims, pos, splits, slice, 'E', os.E, scenario.E_floating, [0,0])
    setBoundsVar(bounds, fullDims, pos, splits, slice, 'U', os.U, scenario.U_floating, [1,0])
    setBoundsVar(bounds, fullDims, pos, splits, slice, 'V', os.V, scenario.V_floating, [0,1])
    setBoundsVar(bounds, fullDims, pos, splits, slice, 'S', os.S, False, [0, 0])
    setBoundsVar(bounds, fullDims, pos, splits, slice, 'T', os.T, False, [0, 0])
    if sp.passiveTracer:
        setBoundsVar(bounds, fullDims, pos, splits, slice, 'X', os.X, False, [0, 0])


# Splits the given bound matrix in a baroclinic (3D) and barotropic (depth-averaged) part.
# Bnd is expected to have two dimensions, with the first being horizontal and the second being vertical.
def splitUVBound(dwd, bnd):

    ua = np.zeros((bnd.shape[0]))
    hua = np.zeros((bnd.shape[0]))
    ub = np.zeros(bnd.shape)
    sumD = np.zeros((bnd.shape[0]))
    for k in range(0, dwd.shape[1]):
        sumD = sumD + dwd[:,k]
        hua = hua + np.multiply(dwd[:,k], bnd[:,k])
    sumD[sumD==0] = 1e-6
    ua = np.divide(hua, sumD)
    for k in range(0, dwd.shape[1]):
        ub[:,k] = bnd[:,k] - ua
    return hua, ua, ub


# Apply the bounds for a given variable retrieved from the scenario (for the full domain)
# to our slice of the domain. This function is used specifically for U or V in a mode split
# setting where the U/V values need to be split into baroclinic/barotropic modes.
def setBoundsSplitVar(bounds, fullDims, pos, splits, slice, varName, varA, varHA, varB, floating, os):
    # First check if bounds have been specified for this variable:
    if not varName in bounds:
        if not floating:
            return
        else:
            # No bounds given, but "floating" is set to true. This means we should
            # copy the values right inside the edges out to the edges:
            copyInnerValues(pos, splits, varA)
            copyInnerValues(pos, splits, varHA)
            copyInnerValues(pos, splits, varB)
            return

    bnd = bounds[varName] # [left, bottom, right, upper]

    # Left:
    if pos[0]==0:
        if varName=='U':
            bHA, bA, bB = splitUVBound(os.DWD[0,:,:], bnd[0][slice[2]:slice[3],...])
            os.UA[0,:] = bA[...]
            os.HUA[0, :] = bHA[...]
            os.UB[0,:,:] = bB[...]

        elif varName=='V':
            bHA, bA, bB = splitUVBound(os.DSD[0, :, :], bnd[0][slice[2]:slice[3]-1, ...])
            os.VA[0, :] = bA[...]
            os.HVA[0, :] = bHA[...]
            os.VB[0, :, :] = bB[...]

    # Bottom:
    if pos[1]==0:
        if varName=='U':
            bHA, bA, bB = splitUVBound(os.DWD[:,0,:], bnd[1][slice[0]:slice[1]-1,...])
            os.UA[:,0] = bA[...]
            os.HUA[:,0] = bHA[...]
            os.UB[:,0,:] = bB[...]

        elif varName=='V':
            bHA, bA, bB = splitUVBound(os.DSD[:,0, :], bnd[1][slice[0]:slice[1], ...])
            os.VA[:,0] = bA[...]
            os.HVA[:,0] = bHA[...]
            os.VB[:,0, :] = bB[...]

    # Right:
    if pos[0]==splits[0]-1:
        if varName == 'U':
            bHA, bA, bB = splitUVBound(os.DWD[-1, :, :], bnd[2][slice[2]:slice[3], ...])
            os.UA[-1, :] = bA[...]
            os.HUA[-1, :] = bHA[...]
            os.UB[-1, :, :] = bB[...]

        elif varName == 'V':
            bHA, bA, bB = splitUVBound(os.DSD[-1, :, :], bnd[2][slice[2]:slice[3]-1, ...])
            os.VA[-1, :] = bA[...]
            os.HVA[-1, :] = bHA[...]
            os.VB[-1, :, :] = bB[...]

    # Top:
    if pos[1]==splits[1]-1:
        if varName=='U':
            bHA, bA, bB = splitUVBound(os.DWD[:,-1,:], bnd[3][slice[0]:slice[1]-1,...])
            os.UA[:,-1] = bA[...]
            os.HUA[:,-1] = bHA[...]
            os.UB[:,-1,:] = bB[...]

        elif varName=='V':
            bHA, bA, bB = splitUVBound(os.DSD[:,-1, :], bnd[3][slice[0]:slice[1], ...])
            os.VA[:,-1] = bA[...]
            os.HVA[:,-1] = bHA[...]
            os.VB[:,-1, :] = bB[...]


# Update boundary values for E, U and V states as part of the mode-splitting loop.
def setEUVBounds(scenario, fullDims, fullDepth, pos, splits, slice, t, os, sp):
    if pos[0]>0 and pos[0]<splits[0]-1 and pos[1]>0 and pos[1]<splits[1]-1:
        return # This slice does not have any external boundaries
    bounds = scenario.setBounds(fullDims[0], fullDims[1], fullDims[2], fullDepth, t, os)

    setBoundsVar(bounds, fullDims, pos, splits, slice, 'E', os.E, scenario.E_floating, [0,0])
    setBoundsSplitVar(bounds, fullDims, pos, splits, slice, 'U', os.UA, os.HUA, os.UB, scenario.U_floating, os)
    setBoundsSplitVar(bounds, fullDims, pos, splits, slice, 'V', os.VA, os.HVA, os.VB, scenario.V_floating, os)
