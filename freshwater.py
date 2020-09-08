

def addFreshwater(scenario, fullDims, pos, splits, slice, t, os, sp):
    riverPos, riverDir, flow = scenario.getFreshwater(fullDims[0], fullDims[1], t, os)
    for riv in range(0, len(riverPos)):
        pos = riverPos[riv]

        if pos[0]>=slice[0] and pos[0]<slice[1] and pos[1]>=slice[2] and pos[1]<slice[3]:
            # River is within our slice. Find its position in our slice:
            i = pos[0]-slice[0]
            j = pos[1]-slice[2]
            rFlow = flow[riv]
            flowSpeed = rFlow[0]/sp.dx/os.cellHeights[i,j,0]
            if riverDir[riv]=="west":
                os.cellHeights[i-1,j,0] = os.cellHeights[i,j,0]
                if sp.modeSplittingOn:
                    os.HUA[i-1,j,0] = flowSpeed*os.cellHeights[i,j,0]
                    os.DW[i-1,j,0] = os.cellHeights[i,j,0]
                else:
                    os.U[i-1,j,0] = flowSpeed
                os.T[i-1,j,0] = rFlow[1]
                os.S[i-1,j,0] = rFlow[2]
            elif riverDir[riv] == "south":
                os.cellHeights[i,j-1,0] = os.cellHeights[i,j,0]
                if sp.modeSplittingOn:
                    os.HVA[i,j-1,0] = flowSpeed*os.cellHeights[i,j,0]
                    os.DS[i,j-1,0] = os.cellHeights[i,j,0]
                else:
                    os.V[i,j-1,0] = flowSpeed
                os.T[i,j-1,0] = rFlow[1]
                os.S[i,j-1,0] = rFlow[2]
            elif riverDir[riv] == "east":
                os.cellHeights[i+1,j,0] = os.cellHeights[i,j,0]
                if sp.modeSplittingOn:
                    os.HUA[i,j] = -flowSpeed*os.cellHeights[i,j,0]
                    os.DW[i,j,0] = os.cellHeights[i,j,0]
                else:
                    os.U[i,j,0] = -flowSpeed
                os.T[i+1,j,0] = rFlow[1]
                os.S[i+1,j,0] = rFlow[2]
            elif riverDir[riv] == "north":
                os.cellHeights[i,j+1,0] = os.cellHeights[i,j,0]
                if sp.modeSplittingOn:
                    os.HVA[i,j,0] = -flowSpeed*os.cellHeights[i,j,0]
                    os.DS[i,j,0] = os.cellHeights[i,j,0]
                else:
                    os.V[i,j,0] = -flowSpeed
                os.T[i,j+1,0] = rFlow[1]
                os.S[i,j+1,0] = rFlow[2]
            else: # top
                print("Wrong river direction specifier: "+str(riverDir[riv]))

    return