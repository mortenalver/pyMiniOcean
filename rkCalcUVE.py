import numpy as np
from oceanState import *
import math
import threading
from utils import getUV, calcMean
from superbeeAdv import superbeeAdv
from biharmon import biharmon
from verticalSpeeds import *


def rkCalcUVE(os, sp, t, dt):
    imax = os.imax
    jmax = os.jmax
    kmax = os.kmax

    f_U = np.zeros(os.U.shape)
    f_V = np.zeros(os.V.shape)
    f_E = np.zeros(os.E.shape)

    # Based on the values in os.W[i,j,0] - the vertical speeds calculated at the top edge
    # of the surface layer cells, calculate the updated elevation:
    computeVerticalSpeeds(os, sp)
    f_E[1:-1,1:-1] = os.W[1:-1,1:-1,0]

    # Create temporary arrays:
    p_above = sp.p_atm0*np.ones((imax,jmax))
    p_diff = np.zeros((imax,jmax))
    p_gradx = np.zeros((imax-1,jmax,kmax))
    p_grady = np.zeros((imax, jmax-1, kmax))

    # Calculate pressure gradients:
    for k in range(0,kmax):
        for i in range(0,imax):
            for j in range(0,jmax):
                if k<os.kmm[i,j]:
                    p_diff[i,j] = os.cellHeights[i,j,k]*9.81*os.rho[i,j,k]
                else:
                    p_diff[i,j] = 0


        for i in range(0,imax-1):
            for j in range(0,jmax):
                if os.maskU[i,j,k]>0:
                    # height measured from below:
                    meanCellHeight = 0.5 * (os.cellHeights[i,j,k] + os.cellHeights[i+1,j,k]);
                    if k==0: # Surface layer
                        p_gradx[i,j,k] = (p_above[i+1,j] + p_diff[i+1,j]*(os.cellHeights[i+1,j,k]-0.5*meanCellHeight)/os.cellHeights[i+1,j,k] \
                             - (p_above[i,j] + p_diff[i,j]*(os.cellHeights[i,j,k]-0.5*meanCellHeight)/os.cellHeights[i,j,k])) / sp.dx;
                    else: # Mid or bottom layer:
                        p_gradx[i,j,k] = (p_above[i+1,j] + p_diff[i+1,j]*0.5*meanCellHeight/os.cellHeights[i+1,j,k] \
                            - (p_above[i,j] + p_diff[i,j]* 0.5*meanCellHeight/os.cellHeights[i,j,k])) / sp.dx




        for i in range(0,imax):
            for j in range(0,jmax-1):
                if os.maskV[i,j,k]>0:
                    # height measured from below:
                    meanCellHeight = 0.5 * (os.cellHeights[i,j,k] + os.cellHeights[i,j+1,k]);
                    if k==0: # Surface layer
                        p_grady[i,j,k] = (p_above[i,j+1] + p_diff[i,j+1]*(os.cellHeights[i,j+1,k]-0.5*meanCellHeight)/os.cellHeights[i,j+1,k] \
                             - (p_above[i,j] + p_diff[i,j]*(os.cellHeights[i,j,k]-0.5*meanCellHeight)/os.cellHeights[i,j,k])) / sp.dx;
                    else: # Mid or bottom layer:
                        p_grady[i,j,k] = (p_above[i,j+1] + p_diff[i,j+1]*0.5*meanCellHeight/os.cellHeights[i,j+1,k] \
                            - (p_above[i,j] + p_diff[i,j]* 0.5*meanCellHeight/os.cellHeights[i,j,k])) / sp.dx

        p_above = p_above + p_diff


    # Calculate horizontal accelerations in U direction:
    for i in range(0,imax-1):
        for j in range(0,jmax):
            if os.kmm[i,j] >= 0:
                for k in range(0, os.kmax):

                    if not os.maskU[i,j,k]:
                        f_U[i,j,k:os.kmax] = math.nan
                        break

                    # Get the value at this point:
                    val = os.U[i,j,k]

                    # Estimate the local V and W values by interpolation:
                    vMean = calcMean((getUV(os.V,i,j-1,k,math.nan), getUV(os.V,i,j,k,math.nan),
                            getUV(os.V,i+1,j-1, k, math.nan), getUV(os.V,i+1,j,k,math.nan)))
                    wMean = calcMean((getUV(os.W,i,j,k, math.nan), getUV(os.W,i+1,j,k, math.nan)))

                    # Estimate the local d2u/dz2 (double derivative):
                    if k>0:
                        dz_up = 0.5*(os.cellHeights[i,j,k]+os.cellHeights[i,j,k-1])
                    else:
                        dz_up = os.cellHeights[i,j,k]
                    if k<kmax-1:
                        dz_down = 0.5*(os.cellHeights[i,j,k]+os.cellHeights[i,j,k+1])
                    else:
                        dz_down = os.cellHeights[i,j,k]
                    d2u_dz2 = ((getUV(os.U,i,j,k-1,val) - val)/dz_up \
                               - (val - getUV(os.U,i,j,k+1,val))/dz_down)/(0.5*(dz_up+dz_down))

                    if sp.biharmonic:
                        # If biharmonic is activated the diffusion is handled later, so we
                        # can set it to 0 for now:
                        diffUV = 0
                    else:
                        # Estimate the local d2u/dx2 (double derivative):
                        d2u_dx2 = (getUV(os.U,i-1,j,k,val) - 2*val + getUV(os.U,i+1,j,k,val))/(sp.dx*sp.dx)
                        # Estimate the local d2u/dy2 (double derivative):
                        d2u_dy2 = (getUV(os.U,i,j-1,k,val) - 2*val + getUV(os.U,i,j+1,k,val))/(sp.dx*sp.dx)
                        # Calculate diffusion term:
                        diffUV = os.AH[i,j,k]*(d2u_dx2 + d2u_dy2)

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



                    # Sum up the terms and calculate next time step value:
                    f_U[i,j,k] = (
                        - p_gradx[i,j,k]/sp.rho_0  # Pressure term
                        + advU + advV + advW   # Advective terms
                        + sp.A_z*d2u_dz2   # Vertical eddy viscosity
                        + diffUV # Horizontal diffusion
                        + 2*sp.omega*math.sin(sp.phi0)*vMean) # Coriolis




    # Calculate horizontal accelerations in V direction:
    for i in range(0,imax):
        for j in range(0,jmax-1):
            if os.kmm[i,j] >= 0:
                for k in range(0, os.kmm[i,j]):

                    if not os.maskV[i,j,k]:
                        f_V[i,j,k:os.kmax] = math.nan
                        break

                    # Get the value at this point:
                    val = os.V[i, j, k]

                    # Estimate the local U and W values by interpolation:
                    uMean = calcMean((getUV(os.U,i-1,j,k,math.nan), getUV(os.U,i,j,k,math.nan),
                            getUV(os.U, i-1, j+1,k,math.nan), getUV(os.U,i,j+1,k,math.nan)))
                    wMean = calcMean((getUV(os.W, i, j, k, math.nan), getUV(os.W, i, j+1, k, math.nan)))

                    # Estimate the local d2u/dz2 (double derivative):
                    if k>0:
                        dz_up = 0.5*(os.cellHeights[i,j,k]+os.cellHeights[i,j,k-1])
                    else:
                        dz_up = os.cellHeights[i,j,k]
                    if k<kmax-1:
                        dz_down = 0.5*(os.cellHeights[i,j,k]+os.cellHeights[i,j,k+1])
                    else:
                        dz_down = os.cellHeights[i,j,k]
                    d2u_dz2 = ((getUV(os.V,i,j,k-1,val) - val)/dz_up \
                               - (val - getUV(os.V,i,j,k+1,val))/dz_down)/(0.5*(dz_up+dz_down))

                    if sp.biharmonic:
                        # If biharmonic is activated the diffusion is handled later, so we
                        # can set it to 0 for now:
                        diffUV = 0
                    else:
                        # Estimate the local d2u/dx2 (double derivative):
                        d2u_dx2 = (getUV(os.V,i-1,j,k,val) - 2*val + getUV(os.V,i+1,j,k,val))/(sp.dx*sp.dx)
                        # Estimate the local d2u/dy2 (double derivative):
                        d2u_dy2 = (getUV(os.V,i,j-1,k,val) - 2*val + getUV(os.V,i,j+1,k,val))/(sp.dx*sp.dx)
                        # Calculate diffusion term:
                        diffUV = os.AH[i,j,k]*(d2u_dx2 + d2u_dy2)

                    # Calculate nonlinear (advective) terms:
                    if sp.advectiveTermsOn:
                        # Calculate the advection (nonlinear) terms using the
                        # Superbee flux limiter to limit oscillations while
                        # suppressing numerical diffusion:
                        advU = superbeeAdv(dt, sp.dx, getUV(os.V,i-2,j,k,val), getUV(os.V,i-1,j,k,val), val,
                            getUV(os.V,i+1,j,k,val), getUV(os.V,i+2,j,k,val), val, val)
                        advV = superbeeAdv(dt, sp.dx, getUV(os.V,i,j-2,k,val), getUV(os.V,i,j-1,k,val), val,
                            getUV(os.V,i,j+1,k,val), getUV(os.V,i,j+2,k,val), vMean, vMean)
                        advW = superbeeAdv(dt, sp.dx, getUV(os.V,i,j,k-2,val), getUV(os.V,i,j,k-1,val), val,
                            getUV(os.V,i,j,k+1,val), getUV(os.V,i,j+2,k+2,val), wMean, wMean)
                    else:
                        advU = 0
                        advV = 0
                        advW = 0

                    # Sum up the terms and calculate next time step value:
                    f_V[i,j,k] = (
                        - p_grady[i,j,k]/sp.rho_0   # Pressure term
                        + advU + advV + advW   # Advective terms
                        + sp.A_z * d2u_dz2   # Vertical eddy viscosity
                        + diffUV  # Horizontal diffusion
                        -2*sp.omega*math.sin(sp.phi0)*uMean)  # Coriolis



    # If we are using biharmonic diffusion of velocities, do it here:
    if sp.biharmonic:
        (diffU, diffV) = biharmon(os, sp, os.U, os.V)
        f_U = f_U - diffU
        f_V = f_V - diffV


    # Wind stress and bottom friction, U:
    for i in range(0,os.imax-1):
        for j in range(0,os.jmax):
            if os.maskU[i,j,0] > 0: # Check if there is a valid current vector at this position
                # Surface cell, average height on cell border:
                dz_mean = 0.5*(os.cellHeights[i,j,0]+os.cellHeights[i+1,j,0])
                f_U[i,j,0] = f_U[i,j,0] + sp.windStressCoeff*os.windU[i,j]/dz_mean

                # Bottom friction. Apply at the minimum kmax of the
                # neighbouring cells. We need to calculate the absolute value
                # of the current speed here, based on U and interpolated V values:
                k = min(os.kmm[i,j], os.kmm[i+1,j])-1
                # Bottom cell, average height on cell border:
                dz_mean = 0.5*(os.cellHeights[i,j,k]+os.cellHeights[i+1,j,k])
                # V value interpolated here:
                meanV = calcMean([getUV(os.V,i,j-1,k,math.nan), getUV(os.V,i+1,j-1,k,math.nan),
                            getUV(os.V,i,j,k,math.nan), getUV(os.V,i+1,j,k,math.nan)])
                speed = math.sqrt(os.U[i,j,k]*os.U[i,j,k] + meanV*meanV);
                f_U[i,j,k] = f_U[i,j,k] - sp.C_b*os.U[i,j,k]*speed/dz_mean

    # Wind stress and bottom friction, V:
    for i in range(0,os.imax):
        for j in range(0,os.jmax-1):
            if os.maskV[i,j,0] > 0: # Check if there is a valid current vector at this position
                # Surface cell, average height on cell border:
                dz_mean = 0.5*(os.cellHeights[i,j,0]+os.cellHeights[i,j+1,0])
                f_V[i, j, 0] = f_V[i, j, 0] + sp.windStressCoeff*os.windV[i, j]/dz_mean

                # Bottom friction. Apply at the minimum kmax of the
                # neighbouring cells. We need to calculate the absolute value
                # of the current speed here, based on U and interpolated V values:
                k = min(os.kmm[i,j], os.kmm[i,j+1])-1
                # Bottom cell, average height on cell border:
                dz_mean = 0.5*(os.cellHeights[i,j,k]+os.cellHeights[i,j+1,k])
                # V value interpolated here:
                meanU = calcMean([getUV(os.U,i-1,j,k,math.nan), getUV(os.U,i-1,j+1,k,math.nan),
                            getUV(os.U,i,j,k,math.nan), getUV(os.U,i,j+1,k,math.nan)])
                speed = math.sqrt(os.V[i,j,k]*os.V[i,j,k] + meanV*meanV);
                f_V[i,j,k] = f_V[i,j,k] - sp.C_b*os.V[i,j,k]*speed/dz_mean



    #for i in range(0,imax-1):
    #    for j in range(0,jmax):
    #        if os.kmm[i,j] >= 0:
    #            for k in range(0, os.kmax):
    #                if os.maskU[i, j, k] and math.isnan(f_U[i,j,k]):
    #                    print((i,j,k))
    #                    print((i, j, k))


    return f_U, f_V, f_E
