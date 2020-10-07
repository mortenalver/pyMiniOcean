import numpy as np
from verticalSpeeds import computeVerticalSpeeds
from utils import getUV, getUV2, calcMean, vertAverageUBVB
from superbeeAdv import superbeeAdv
from biharmon import biharmon, biharmon2D
import setBounds
import math
import mpi
import matplotlib.pyplot as plt

# Mode splitting integration scheme.
def integrate(os, sp, scenario, fullDims, fullDepth, pos, splits, slice, t, doMpi, comm, rank):
    # Some small precalculations:
    dt = sp.dt/sp.nsub
    dtn = sp.dt
    dx = sp.dx
    dx2 = dx*dx
    dtdx = dt/dx
    dtndx = dtn/dx

    setBounds.setEUVBounds(scenario, fullDims, fullDepth, pos, splits, slice, t, os, sp)

    UB_prelim = os.UB.copy()
    VB_prelim = os.VB.copy()
    AU = np.zeros((os.imax-1, os.jmax))
    AV = np.zeros((os.imax, os.jmax-1))

    # If we are using biharmonic diffusion of velocities, compute those for velocity deviations here:
    if sp.biharmonic:
        (diffU, diffV) = biharmon(os, sp, os.UB, os.VB)


    # Then calculate preliminary U deviations for next long time step by including
    # advection (of full velocities), coriolis terms based on deviations, pressure term based on
    # deviations and diffusion of deviations:
    for i in range(0, os.imax-1):
        for j in range(0, os.jmax):
            # First define the number of layers here:
            kmx = min(os.kmm[i,j], os.kmm[i+1,j])
            sumD = 0
            for k in range(0, kmx):
                #if not os.maskU[i,j,k]:
                #    break
                # As we go down through the layers we sum up the pressure/rho values in the neighbouring cells:
                if k==0:
                    # Center of cell boundary:
                    p_C = np.array([9.81*0.5*os.DW[i,j,k]*(os.rho[i,j,k]/sp.rho_0 - 1),
                          9.81*0.5*os.DW[i,j,k]*(os.rho[i+1,j,k]/sp.rho_0 - 1)])
                    # Bottom of cell boundary, stored for next iteration:
                    p_B = 2*p_C
                else:
                    # Added pressure/rho in this layer:
                    p_h = np.array([9.81*os.DW[i,j,k]*(os.rho[i,j,k]/sp.rho_0 - 1),
                           9.81*os.DW[i,j,k]*(os.rho[i+1,j,k]/sp.rho_0 - 1)])
                    # Center of cell boundary:
                    p_C = p_B + 0.5*p_h
                    # Bottom of cell boundary, stored for next iteration:
                    p_B = p_B + p_h


                # Add pressure term:
                UB_prelim[i,j,k] = UB_prelim[i,j,k] - dtndx*(p_C[1]-p_C[0])

                #if i==10 and j==10:
                #    print("pdiff="+str(p_C[1]-p_C[0])+", p_C[0]="+str(p_C[0])+", sumD="+str(sumD+os.DWD[i,j,k])+
                #          ", rho="+str(os.rho[i,j,k])+", rho/rho0-1="+str(os.rho[i,j,k]/sp.rho_0 - 1))
                #    if k==os.kmax-1:
                #        print()

                # Get the full speed and deviation at this point:
                val = os.UA[i,j] + os.UB[i,j,k]
                valB = os.UB[i,j,k]
                # Estimate the local full V and W values by interpolation:
                vMean = calcMean((getUV(os.V,i,j-1,k,math.nan), getUV(os.V,i,j,k,math.nan),
                        getUV(os.V,i+1,j-1, k, math.nan), getUV(os.V,i+1,j,k,math.nan)))
                wMean = calcMean((getUV(os.W,i,j,k, math.nan), getUV(os.W,i+1,j,k, math.nan)))
                # Estimate the local V deviation value by interpolation:
                vbMean = calcMean((getUV(os.VB,i,j-1,k,math.nan), getUV(os.VB,i,j,k,math.nan),
                        getUV(os.VB,i+1,j-1, k, math.nan), getUV(os.VB,i+1,j,k,math.nan)))
                # Get the absolute speed value:
                absSpeed = math.sqrt(val*val + vMean*vMean)

                # Calculate nonlinear (advective) terms:
                if sp.advectiveTermsOn:
                    # Calculate the advection (nonlinear) terms using the
                    # Superbee flux limiter to limit oscillations while
                    # suppressing numerical diffusion:
                    advU = superbeeAdv(dt, sp.dx, getUV(os.U,i-2,j,k,val), getUV(os.U,i-1,j,k,val), val,
                        getUV(os.U,i+1,j,k,val), getUV(os.U,i+2,j,k,val), val, val)
                    advV = superbeeAdv(dt, sp.dx, getUV(os.U,i,j-2,k,val), getUV(os.U,i,j-1,k,val), val,
                        getUV(os.U,i,j+1,k,val), getUV(os.U,i,j+2,k,val), vMean, vMean)
                    advW = superbeeAdv(dt, sp.dx, getUV(os.U,i,j,k-2,val), getUV(os.U,i,j,k-1,val), val,
                        getUV(os.U,i,j,k+1,val), getUV(os.U,i,j+2,k+2,val), wMean, wMean)
                else:
                    advU = 0
                    advV = 0
                    advW = 0

                # Add advective terms:
                UB_prelim[i,j,k] = UB_prelim[i,j,k] + dtn*(advU + advV + advW)

                # Diffusion:
                if sp.biharmonic:
                    UB_prelim[i,j,k] = UB_prelim[i,j,k] - dtn*diffU[i,j,k]
                else:
                    # Estimate the local d2u/dx2 (double derivative):
                    d2u_dx2 = (getUV(os.UB,i-1,j,k,valB) - 2*valB + getUV(os.UB,i+1,j,k,valB))/dx2
                    # Estimate the local d2u/dy2 (double derivative):
                    d2u_dy2 = (getUV(os.UB,i,j-1,k,valB) - 2*valB + getUV(os.UB,i,j+1,k,valB))/dx2
                    # Calculate diffusion term:
                    UB_prelim[i, j, k] = UB_prelim[i, j, k] + dtn*os.AH[i,j,k]*(d2u_dx2 + d2u_dy2)

                # Coriolis:
                UB_prelim[i,j,k] = UB_prelim[i,j,k] + dtn*2*sp.omega*math.sin(sp.phi0)*vbMean

                # Wind stress. Only if this is the surface layer:
                if k==0:
                    UB_prelim[i,j,k] = UB_prelim[i,j,k] + dtn*sp.windStressCoeff*os.windU[i, j]/os.DWD[i,j,k]

                # Bottom friction. Only if this is bottom layer:
                if k==kmx-1:
                    UB_prelim[i,j,k] = UB_prelim[i,j,k] - dtn*sp.C_b*val*absSpeed/os.DWD[i,j,k]

                # Vertical diffusion.
                if kmx>1 and not sp.implicitVerticalDiff:
                    # Estimate the local d2u/dz2 (double derivative):
                    if k>0:
                        dz_up = 0.5*(os.DWD[i,j,k]+os.DWD[i,j,k-1])
                    else:
                        dz_up = os.DWD[i,j,k]
                    if k<os.kmax-1 and os.maskU[i,j,k+1]:
                        dz_down = 0.5*(os.DWD[i,j,k]+os.DWD[i,j,k+1])
                    else:
                        dz_down = os.DW[i,j,k]
                    d2u_dz2 = ((getUV(os.UB,i,j,k-1,valB) - valB)/dz_up \
                               - (valB - getUV(os.UB,i,j,k+1,valB))/dz_down)/(0.5*(dz_up+dz_down))
                    UB_prelim[i, j, k] = UB_prelim[i, j, k] + dtn * sp.A_z * d2u_dz2

                if math.isnan(UB_prelim[i, j, k]):
                    print(i)


                # Find AU (depth integrated UB for this time step):
                AU[i,j] = AU[i,j] + os.DWD[i,j,k]*UB_prelim[i,j,k]
                sumD = sumD + os.DWD[i,j,k]


            # Divide AU by sum depth and time step:
            if sumD > 0:
                AU[i,j] = AU[i,j]/(dtn*sumD)

                # Subtract AU in all layers:
                for k in range(0, kmx):
                    UB_prelim[i,j,k] = UB_prelim[i,j,k] - AU[i,j]*dtn

            # Vertical diffusion - implicit calculation:
            if sp.implicitVerticalDiff and kmx>1:
                AP = np.zeros((kmx,))
                CP = np.zeros((kmx,))
                SP = np.zeros((kmx,))
                EP = np.zeros((kmx,))
                AP[0] = 0
                for k in range(1,kmx):
                    AP[k] = sp.A_z*dtn/(os.DWD[i,j,k]*0.5*(os.DWD[i,j,k]+os.DWD[i,j,k-1]))
                    CP[k-1] = sp.A_z*dtn/(os.DWD[i,j,k-1]*0.5*(os.DWD[i,j,k-1]+os.DWD[i,j,k]))
                CP[kmx-1] = 0

                SP[0] = 1 + CP[0] + AP[0]
                EP[0] = UB_prelim[i,j,0]
                for k in range(1, kmx):
                    SP[k] = 1 + AP[k] + CP[k] - AP[k]*CP[k-1]/SP[k-1]
                    EP[k] = UB_prelim[i,j,k] + AP[k]*EP[k-1]/SP[k-1]
                UB_prelim[i,j,kmx-1] = EP[kmx-1]/SP[kmx-1]
                for k in range(kmx-2, -1, -1):
                    UB_prelim[i,j,k] = (EP[k]+CP[k]*UB_prelim[i,j,k+1])/SP[k]


    # Calculate preliminary V deviations;
    for i in range(0, os.imax):
        for j in range(0, os.jmax-1):
            # First define the number of layers here:
            kmx = min(os.kmm[i,j], os.kmm[i,j+1])
            sumD = 0
            for k in range(0, kmx):
                #if not os.maskV[i,j,k]:
                #    break
                # As we go down through the layers we sum up the pressure/rho values in the neighbouring cells:
                if k==0:
                    # Center of cell boundary:
                    p_C = np.array([9.81*0.5*os.DS[i,j,k]*(os.rho[i,j,k]/sp.rho_0 - 1),
                          9.81*0.5*os.DS[i,j,k]*(os.rho[i,j+1,k]/sp.rho_0 - 1)])
                    # Bottom of cell boundary, stored for next iteration:
                    p_B = 2*p_C
                else:
                    # Added pressure/rho in this layer:
                    p_h = np.array([9.81*os.DS[i,j,k]*(os.rho[i,j,k]/sp.rho_0 - 1),
                           9.81*os.DS[i,j,k]*(os.rho[i,j+1,k]/sp.rho_0 - 1)])
                    # Center of cell boundary:
                    p_C = p_B + 0.5*p_h
                    # Bottom of cell boundary, stored for next iteration:
                    p_B = p_B + p_h



                # Add pressure term:
                VB_prelim[i,j,k] = VB_prelim[i,j,k] - dtndx*(p_C[1]-p_C[0])

                # Get the full speed and deviation at this point:
                val = os.V[i,j,k]
                valB = os.VB[i,j,k]
                # Estimate the local U and W values by interpolation:
                uMean = calcMean((getUV(os.U,i-1,j,k,math.nan), getUV(os.U,i,j,k,math.nan),
                        getUV(os.U, i-1, j+1,k,math.nan), getUV(os.U,i,j+1,k,math.nan)))
                wMean = calcMean((getUV(os.W, i, j, k, math.nan), getUV(os.W, i, j+1, k, math.nan)))
                # Estimate the local U deviation value by interpolation:
                ubMean = calcMean((getUV(os.UB,i-1,j,k,math.nan), getUV(os.UB,i,j,k,math.nan),
                        getUV(os.UB, i-1, j+1,k,math.nan), getUV(os.UB,i,j+1,k,math.nan)))
                # Get the absolute speed value:
                absSpeed = math.sqrt(val*val + uMean*uMean)

                # Calculate nonlinear (advective) terms:
                if sp.advectiveTermsOn:
                    # Calculate the advection (nonlinear) terms using the
                    # Superbee flux limiter to limit oscillations while
                    # suppressing numerical diffusion:
                    advU = superbeeAdv(dt, sp.dx, getUV(os.V,i-2,j,k,val), getUV(os.V,i-1,j,k,val), val,
                        getUV(os.V,i+1,j,k,val), getUV(os.V,i+2,j,k,val), uMean, uMean)
                    advV = superbeeAdv(dt, sp.dx, getUV(os.V,i,j-2,k,val), getUV(os.V,i,j-1,k,val), val,
                        getUV(os.V,i,j+1,k,val), getUV(os.V,i,j+2,k,val), val, val)
                    advW = superbeeAdv(dt, sp.dx, getUV(os.V,i,j,k-2,val), getUV(os.V,i,j,k-1,val), val,
                        getUV(os.V,i,j,k+1,val), getUV(os.V,i,j+2,k+2,val), wMean, wMean)
                else:
                    advU = 0
                    advV = 0
                    advW = 0

                # Add advective terms:
                VB_prelim[i,j,k] = VB_prelim[i,j,k] + dtn*(advU + advV + advW)

                # Diffusion:
                if sp.biharmonic:
                    VB_prelim[i,j,k] = VB_prelim[i,j,k] - dtn*diffV[i,j,k]
                else:
                    # Estimate the local d2v/dx2 (double derivative):
                    d2v_dx2 = (getUV(os.VB,i-1,j,k,valB) - 2*valB + getUV(os.VB,i+1,j,k,valB))/dx2
                    # Estimate the local d2v/dy2 (double derivative):
                    d2v_dy2 = (getUV(os.VB,i,j-1,k,valB) - 2*valB + getUV(os.VB,i,j+1,k,valB))/dx2
                    # Calculate diffusion term:
                    VB_prelim[i, j, k] = VB_prelim[i, j, k] + dtn*os.AH[i,j,k]*(d2v_dx2 + d2v_dy2)

                # Coriolis:
                VB_prelim[i,j,k] = VB_prelim[i,j,k] - dtn*2*sp.omega*math.sin(sp.phi0)*ubMean

                # Wind stress. Only if this is the surface layer:
                if k==0:
                    VB_prelim[i,j,k] = VB_prelim[i,j,k] + dtn*sp.windStressCoeff*os.windV[i,j]/os.DSD[i,j,k]

                # Bottom friction. Only if this is bottom layer:
                if k == kmx-1:
                    VB_prelim[i, j, k] = VB_prelim[i,j,k] - sp.C_b*val*absSpeed/os.DSD[i,j,k]

                # Vertical diffusion.
                if kmx>1 and not sp.implicitVerticalDiff:
                    # Estimate the local d2u/dz2 (double derivative):
                    if k > 0:
                        dz_up = 0.5 * (os.DSD[i,j,k] + os.DSD[i,j,k-1])
                    else:
                        dz_up = os.DSD[i,j,k]
                    if k<os.kmax-1 and os.maskV[i,j,k+1]:
                        dz_down = 0.5 * (os.DSD[i,j,k] + os.DSD[i,j,k+1])
                    else:
                        dz_down = os.DS[i,j,k]
                    d2u_dz2 = ((getUV(os.VB,i,j,k-1,valB)-valB)/dz_up
                               - (valB-getUV(os.VB,i,j,k+1, valB)) / dz_down) / (0.5 * (dz_up + dz_down))
                    VB_prelim[i,j,k] = VB_prelim[i,j,k] + sp.A_z * d2u_dz2

                # Find AV (depth integrated VB for this time step):
                AV[i,j] = AV[i,j] + os.DSD[i,j,k]*VB_prelim[i,j,k]
                sumD = sumD + os.DSD[i,j,k]


            # Divide AV by sum depth and time step:
            if sumD > 0:
                AV[i,j] = AV[i,j]/(dtn*sumD)

                # Subtract AU in all layers:
                for k in range(0, kmx):
                    VB_prelim[i,j,k] = VB_prelim[i,j,k]-AV[i,j]*dtn

            # Vertical diffusion - implicit calculation:
            if sp.implicitVerticalDiff and kmx > 1:
                AP = np.zeros((kmx,))
                CP = np.zeros((kmx,))
                SP = np.zeros((kmx,))
                EP = np.zeros((kmx,))
                AP[0] = 0
                for k in range(1, kmx):
                    AP[k] = sp.A_z*dtn/(os.DSD[i,j,k]*0.5*(os.DSD[i,j,k]+os.DSD[i,j,k-1]))
                    CP[k-1] = sp.A_z*dtn/(os.DSD[i,j,k-1]*0.5*(os.DSD[i,j,k-1]+os.DSD[i,j,k]))
                CP[kmx-1] = 0

                SP[0] = 1 + CP[0] + AP[0]
                EP[0] = VB_prelim[i,j,0]
                for k in range(1, kmx):
                    SP[k] = 1 + AP[k] + CP[k] - AP[k] * CP[k-1] / SP[k-1]
                    EP[k] = VB_prelim[i,j,k] + AP[k] * EP[k-1] / SP[k-1]
                VB_prelim[i,j,kmx-1] = EP[kmx-1] / SP[kmx-1]
                for k in range(kmx-2,-1,-1):
                    VB_prelim[i,j,k] = (EP[k] + CP[k] * VB_prelim[i,j,k+1]) / SP[k]


    # Update 3D speeds:
    os.UB[...] = UB_prelim[...]
    os.VB[...] = VB_prelim[...]

    # Communicate UB/VB between processes:
    #if doMpi:
    #    mpi.communicate3D(comm, rank, pos, splits, os, sp)

    # Vertically integrated model:
    for subt in range(0, sp.nsub):

        # Communicate 2D fields between processes:
        if doMpi:
            mpi.communicate2D(comm, rank, pos, splits, os, sp)

        deltU = np.zeros(os.HUA.shape)
        deltV = np.zeros(os.HVA.shape)

        if sp.biharmonic:
            (diffU, diffV) = biharmon2D(os, sp, os.UA, os.VA)

        # U direction:
        for i in range(0, os.imax-1):
            for j in range(1, os.jmax-1):
                if not os.maskU[i,j,0]:
                    continue

                # Local values:
                kmx = min(os.kmm[i,j], os.kmm[i+1,j])-1 # Number of wet layers
                dwad = os.DWA[i,j] + 0.5*(os.E[i,j] + os.E[i+1,j])
                val = os.UA[i,j]
                hval = os.HUA[i,j]
                # Weighted average of V:
                averHV = 0
                if os.maskV[i,j-1,0]:
                    averHV = averHV + getUV2(os.HVA,i,j-1,0)/os.hvSqr[i,j-1]
                if os.maskV[i+1,j-1,0]:
                    averHV = averHV + getUV2(os.HVA,i+1,j-1,0)/os.hvSqr[i+1,j-1]
                if os.maskV[i,j,0]:
                    averHV = averHV + getUV2(os.HVA,i,j,0)/os.hvSqr[i,j]
                if os.maskV[i+1,j,0]:
                    averHV = averHV + getUV2(os.HVA,i+1,j,0)/os.hvSqr[i+1,j]
                averHV = 0.25*averHV*os.huSqr[i,j]
                # Unweighted:
                averVA = 0.25*(getUV2(os.VA,i,j-1,0)+ getUV2(os.VA,i+1,j-1,0)
                              +getUV2(os.VA,i,j,0)+getUV2(os.VA,i+1,j,0))
                averVB = 0.25*(getUV(os.VB,i,j-1,kmx,0)+ getUV(os.VB,i+1,j-1,kmx,0)
                              +getUV(os.VB,i,j,kmx,0)+getUV(os.VB,i+1,j,kmx,0))
                # Get bottom current to calculate bottom drag:
                ubot = val + os.UB[i,j,kmx]
                vbot = averVA + averVB
                vvec = math.sqrt(ubot*ubot + vbot*vbot)

                # TODO: depth-integrated equation, support for variable atmospheric pressure

                deltU[i,j] = (
                    dt*2*sp.omega*math.sin(sp.phi0)*averHV # Coriolis
                    + dtdx*9.81*dwad*(os.E[i,j]-os.E[i+1,j]) # External gravity waves
                    - dt*sp.C_b*ubot*vvec # Bottom friction
                    + dt*dwad*AU[i,j]) # Average rate of change extracted from 3D step

                if sp.biharmonic:
                    deltU[i,j] = deltU[i,j] - dt*dwad*diffU[i,j]
                else:
                    # Estimate the local d2u/dx2 (double derivative):
                    d2u_dx2 = (getUV2(os.UA,i-1,j,val) - 2*val + getUV2(os.UA,i+1,j,val))/dx2
                    # Estimate the local d2u/dy2 (double derivative):
                    d2u_dy2 = (getUV2(os.UA,i,j-1,val) - 2*val + getUV2(os.UA,i,j+1,val))/dx2
                    # Calculate diffusion term:
                    deltU[i,j] = deltU[i,j] + dt*os.AM2D[i,j]*(d2u_dx2 + d2u_dy2)







        # V direction:
        for i in range(1, os.imax-1):
            for j in range(0, os.jmax-1):
                if not os.maskV[i,j,0]:
                    continue

                # Local values:
                kmx = min(os.kmm[i,j], os.kmm[i,j+1])-1 # Number of wet layers
                dsad = os.DSA[i,j] + 0.5*(os.E[i,j] + os.E[i,j+1])
                val = os.VA[i,j]
                hval = os.HVA[i,j]
                # Weighted average of U:
                averHU = 0
                if os.maskU[i-1,j,0]:
                    averHU = averHU + getUV2(os.HUA,i-1,j,0)/os.huSqr[i-1,j]
                if os.maskU[i,j,0]:
                    averHU = averHU + getUV2(os.HUA,i,j,0)/os.huSqr[i,j]
                if os.maskU[i-1,j+1,0]:
                    averHU = averHU + getUV2(os.HUA,i-1,j+1,0)/os.huSqr[i-1,j+1]
                if os.maskU[i,j+1,0]:
                    averHU = averHU + getUV2(os.HUA,i,j+1,0)/os.huSqr[i,j+1]
                averHU = 0.25*averHU*os.hvSqr[i,j]
                # Unweighted:
                averUA = 0.25*(getUV2(os.UA,i-1,j,0)+ getUV2(os.UA,i,j,0)
                              +getUV2(os.UA,i-1,j+1,0)+getUV2(os.UA,i,j+1,0))
                averUB = 0.25*(getUV(os.UB,i-1,j,kmx,0)+ getUV(os.UB,i,j,kmx,0)
                              +getUV(os.UB,i-1,j+1,kmx,0)+getUV(os.UB,i,j+1,kmx,0))
                # Get bottom current to calculate bottom drag:
                vbot = val + os.VB[i,j,kmx]
                ubot = averUA + averUB
                vvec = math.sqrt(ubot*ubot + vbot*vbot)

                # TODO: depth-integrated equation, support for variable atmospheric pressure

                deltV[i,j] = (
                    - dt*2*sp.omega*math.sin(sp.phi0)*averHU # Coriolis
                    + dtdx*9.81*dsad*(os.E[i,j]-os.E[i,j+1]) # External gravity waves
                    - dt*sp.C_b*vbot*vvec # Bottom friction
                    + dt*dsad*AV[i,j]) # Average rate of change extracted from 3D step

                if sp.biharmonic:
                    deltV[i,j] = deltV[i,j] - dt*dsad*diffV[i,j]
                else:
                    # Estimate the local d2u/dx2 (double derivative):
                    d2u_dx2 = (getUV2(os.VA,i-1,j,val) - 2*val + getUV2(os.VA,i+1,j,val))/dx2
                    # Estimate the local d2u/dy2 (double derivative):
                    d2u_dy2 = (getUV2(os.VA,i,j-1,val) - 2*val + getUV2(os.VA,i,j+1,val))/dx2
                    # Calculate diffusion term:
                    deltV[i,j] = deltV[i,j] + dt*os.AM2D[i,j]*(d2u_dx2 + d2u_dy2)

        # Update HUA:
        os.HUA = os.HUA + deltU

        # Calculate new UA:
        for i in range(0, os.imax - 1):
            for j in range(0, os.jmax):
                if not os.maskU[i, j, 0]:
                    continue
                dwad = os.DWA[i, j] + 0.5 * (os.E[i, j] + os.E[i + 1, j])
                os.UA[i, j] = os.HUA[i, j] / dwad

        # Update HVA:
        os.HVA = os.HVA + deltV

        # Calculate new VA:
        for i in range(0, os.imax):
            for j in range(0, os.jmax-1):
                if not os.maskV[i,j,0]:
                    continue
                dsad = os.DSA[i,j] + 0.5*(os.E[i,j] + os.E[i,j+1])
                os.VA[i, j] = os.HVA[i,j]/dsad

        # Done with short time step for UA/HUA and VA/HVA.
        # Calculate new elevation:
        for i in range(1, os.imax-1):
            for j in range(1, os.jmax-1):
                if os.kmm[i,j] == 0:
                    continue
                os.E[i,j] = os.E[i,j] + dtdx*(os.HUA[i-1,j] - os.HUA[i,j] + os.HVA[i,j-1] - os.HVA[i,j])



        # Update local time variable:
        t = t + dt
        # Update EUV bounds:
        setBounds.setEUVBounds(scenario, fullDims, fullDepth, pos, splits, slice, t, os, sp)


    # Done with depth integrated (short) time steps

    # Recompute U:
    for i in range(0, os.imax-1):
        for j in range(0, os.jmax):
            for k in range(0, os.kmax):
                if not os.maskU[i,j,k]:
                    # Check we are at the surface layer and there is still current. If so, it's a river outlet.
                    if k==0 and os.DW[i,j,0] > 0 and os.HUA[i,j] != 0:
                        os.U[i,j,k] = os.HUA[i,j]/os.DW[i,j,0]
                    break
                os.U[i,j,k] = os.UA[i,j] + os.UB[i,j,k]

    # Recompute V:
    for i in range(0, os.imax):
        for j in range(0, os.jmax-1):
            for k in range(0, os.kmax):
                if not os.maskV[i,j,k]:
                    # Check we are at the surface layer and there is still current. If so, it's a river outlet.
                    if k==0 and os.DS[i,j,0] > 0 and os.HVA[i,j] != 0:
                        os.V[i,j,k] = os.HVA[i,j]/os.DS[i,j,0]
                    break
                os.V[i,j,k] = os.VA[i,j] + os.VB[i,j,k]





    # Recompute vertical speeds:
    computeVerticalSpeeds(os, sp)
