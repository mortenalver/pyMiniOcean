from utils import getUV
import math
import numpy as np

# Advect/diffuse a tracer field based on current speeds, vertical mixing coefficients
# and horizontal mixing coefficients.
def advectField(trc, trc_next, os, sp):

    for i in range(1,os.imax-1):
        for j in range(1,os.jmax-1):
            for k in range(0,os.kmm[i,j]):
                advS = 0
                trc_ij = trc[i,j,k]

                # Find the cell heights interpolated to the cell boundaries:
                vsize = [0.5*(os.cellHeights[i-1,j,k]+os.cellHeights[i,j,k]),
                    0.5*(os.cellHeights[i,j,k]+os.cellHeights[i+1,j,k]),
                    0.5*(os.cellHeights[i,j-1,k]+os.cellHeights[i,j,k]),
                    0.5*(os.cellHeights[i,j,k]+os.cellHeights[i,j+1,k])]

                # For each horizontal direction, if current is into this cell,
                # add advection term of the unit (vol/s)*deltaS:
                if os.maskU[i-1,j,k] and os.U[i-1,j,k] > 0:
                    advS = advS + os.U[i-1,j,k]*sp.dx*vsize[0]*(getUV(trc,i-1,j,k,trc_ij)-getUV(trc,i,j,k,trc_ij))
                if os.maskU[i,j,k] and os.U[i,j,k] < 0:
                    advS = advS - os.U[i,j,k]*sp.dx*vsize[1]*(getUV(trc,i+1,j,k,trc_ij)-getUV(trc,i,j,k,trc_ij))
                if os.maskV[i,j-1,k] and os.V[i,j-1,k] > 0:
                    advS = advS + os.V[i,j-1,k]*sp.dx*vsize[2]*(getUV(trc,i,j-1,k,trc_ij)-getUV(trc,i,j,k,trc_ij))
                if os.maskV[i,j,k] and os.V[i,j,k] < 0:
                    advS = advS - os.V[i,j,k]*sp.dx*vsize[3]*(getUV(trc,i,j+1,k,trc_ij)-getUV(trc,i,j,k,trc_ij))

                # Vertically:
                if os.W[i,j,k+1] > 0:
                    advS = advS + os.W[i,j,k+1]*sp.dx*sp.dx*(getUV(trc,i,j,k+1,trc_ij)-getUV(trc,i,j,k,trc_ij))
                if k > 0:
                    if os.W[i,j,k-1] < 0:
                        advS = advS - os.W[i,j,k-1]*sp.dx*sp.dx*(trc[i,j,k-1]-trc[i,j,k])

                if sp.trcHorizMix:
                    validNb = [0,0,0,0]
                    if os.maskU[i-1,j,k]:
                        validNb[0] = 1
                    if os.maskU[i,j,k]:
                        validNb[1] = 1
                    if os.maskV[i,j-1,k]:
                        validNb[2] = 1
                    if os.maskV[i,j,k]:
                        validNb[3] = 1
                    diffU = (validNb[1]*0.5*(os.AH[i,j,k]+os.AH[i+1,j,k]) *
                        (getUV(trc,i+1,j,k,trc_ij)-getUV(trc,i,j,k,trc_ij))
                        - validNb[0]*0.5*(os.AH[i-1,j,k]+os.AH[i,j,k]) *
                        (getUV(trc,i,j,k,trc_ij)-getUV(trc,i-1,j,k,trc_ij))) / (sp.dx*sp.dx)
                    diffV = (validNb[3]*0.5*(os.AH[i,j,k]+os.AH[i,j+1,k]) *
                        (getUV(trc,i,j+1,k,trc_ij)-getUV(trc,i,j,k,trc_ij))
                        - validNb[2]*0.5*(os.AH[i,j-1,k]+os.AH[i,j,k]) *
                        (getUV(trc,i,j,k,trc_ij)-getUV(trc,i,j-1,k,trc_ij))) / (sp.dx*sp.dx)
                else:
                    diffU = 0
                    diffV = 0

                # # Vertical mixing:
                # if sp.trcVertMix:
                #
                #     if k==0:
                #         v_above = trc_ij
                #         dz_up = os.cellHeights[i,j,k]
                #         kv_above = 0
                #     else:
                #         v_above = trc[i,j,k-1]
                #         dz_up = 0.5*(os.cellHeights[i,j,k]+os.cellHeights[i,j,k-1])
                #         kv_above = os.K_v[i,j,k-1]
                #     if k==os.kmm[i,j]-1:
                #         v_below = trc_ij
                #         dz_down = os.cellHeights[i,j,k]
                #         kv_below = 0
                #     else:
                #         v_below = trc[i,j,k+1]
                #         dz_down = 0.5*(os.cellHeights[i,j,k]+os.cellHeights[i,j,k+1])
                #         kv_below = os.K_v[i,j,k]
                #     diffS = (kv_above*(v_above-trc_ij)/dz_up - kv_below*(trc_ij - v_below)/dz_down)/(0.5*(dz_up+dz_down))
                # else:
                #     diffS = 0

                trc_next[i,j,k] = trc[i,j,k] + sp.dt*(
                            advS/(sp.dx*sp.dx*os.cellHeights[i,j,k]) + diffU + diffV)

                # if i==1 and j==54 and k==5:
                #     print(str(trc_ij)+" -> "+str(trc_next[i,j,k]))
                #if np.isnan(trc_next[i,j,k]):
                #    print("hei")

    # Implicit calculation of vertical mixing of tracers:
    if sp.trcVertMix:
        for i in range(1, os.imax - 1):
            for j in range(1, os.jmax - 1):
                # Vertical diffusion - implicit calculation:
                if os.kmm[i,j] > 1:
                    kmx = os.kmm[i,j]
                    AP = np.zeros((kmx,))
                    CP = np.zeros((kmx,))
                    SP = np.zeros((kmx,))
                    EP = np.zeros((kmx,))
                    AP[0] = 0
                    for k in range(1, kmx):
                        AP[k] = os.K_v[i,j,k-1] * sp.dt/(os.cellHeights[i,j,k]*0.5*(os.cellHeights[i,j,k]+os.cellHeights[i,j,k-1]))
                        CP[k-1] = os.K_v[i,j,k-1]*sp.dt/\
                                  (os.cellHeights[i,j,k-1]*0.5*(os.cellHeights[i,j,k-1]+os.cellHeights[i,j,k]))
                    CP[kmx-1] = 0

                    SP[0] = 1 + CP[0] + AP[0]
                    EP[0] = trc_next[i,j,0]
                    for k in range(1, kmx):
                        SP[k] = 1 + AP[k] + CP[k] - AP[k]*CP[k-1]/SP[k-1]
                        EP[k] = trc_next[i,j,k] + AP[k]*EP[k-1]/SP[k-1]
                    trc_next[i,j,kmx-1] = EP[kmx-1]/SP[kmx-1]
                    for k in range(kmx-2,-1,-1):
                        trc_next[i,j,k] = (EP[k]+CP[k]*trc_next[i,j,k+1])/SP[k]


def advectTempSalt(os, sp):
    advectField(os.T, os.T_next, os, sp)
    advectField(os.S, os.S_next, os, sp)

def advectPassiveTracer(os, sp):
    advectField(os.X, os.X_next, os, sp)
    return