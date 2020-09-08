import math
import numpy as np
from utils import getUV, getUV2

# Calculates biharmonic "diffusion" of velocities
# SINMOD (Sundfjord et al. 2008): "The model calculates horizontal diffusivity of momentum using
# biharmonic friction while horizontal diffusion of scalars applies diffusion
# coefficients as in Smagorinsky (1963)."
def biharmon(os, sp, varU, varV):

    brd=2
    dx4 = sp.dx*sp.dx*sp.dx*sp.dx

    diffU = np.zeros(varU.shape)
    diffV = np.zeros(varV.shape)

    # Diffusion along x direction:
    for i in range(brd, os.imax-1-brd):
        for j in range(brd, os.jmax-brd):
            for k in range(0, os.kmax):
                if not os.maskU[i,j,k]:
                    break
                val = varU[i,j,k]
                s1 = getUV(varU,i+1,j,k,val) + getUV(varU,i,j+1,k,val) \
                    + getUV(varU,i-1,j,k,val) + getUV(varU,i,j-1,k,val)
                s2 = getUV(varU,i+1,j+1,k,val) + getUV(varU,i-1,j+1,k,val) \
                    + getUV(varU,i-1,j-1,k,val) + getUV(varU,i+1,j-1,k,val)
                s3 = os.maskU[i+1,j,k]*getUV(varU,i+2,j,k,val) \
                    + os.maskU[i,j+1,k]*getUV(varU,i,j+2,k,val) \
                    + os.maskU[i-1,j,k]*getUV(varU,i-2,j,k,val) \
                    + os.maskU[i,j-1,k]*getUV(varU,i,j-2,k,val)

                diffU[i,j,k] = sp.KBi*(20*val - 8*s1 + 2*s2 + s3)/dx4

    # Diffusion along y direction:
    for i in range(brd, os.imax-brd):
        for j in range(brd, os.jmax-1-brd):
            for k in range(0, os.kmax):
                if not os.maskV[i,j,k]:
                    break
                val = varV[i,j,k]
                s1 = getUV(varV,i+1,j,k,val) + getUV(varV,i,j+1,k,val) \
                    + getUV(varV,i-1,j,k,val) + getUV(varV,i,j-1,k,val)
                s2 = getUV(varV,i+1,j+1,k,val) + getUV(varV,i-1,j+1,k,val) \
                    + getUV(varV,i-1,j-1,k,val) + getUV(varV,i+1,j-1,k,val)
                s3 = os.maskV[i+1,j,k]*getUV(varV,i+2,j,k,val) \
                    + os.maskV[i,j+1,k]*getUV(varV,i,j+2,k,val) \
                    + os.maskV[i-1,j,k]*getUV(varV,i-2,j,k,val) \
                    + os.maskV[i,j-1,k]*getUV(varV,i,j-2,k,val)

                diffV[i,j,k] = sp.KBi*(20*val - 8*s1 + 2*s2 + s3)/dx4

    return (diffU, diffV)


# Calculates biharmonic "diffusion" of velocities for the depth integrated states
def biharmon2D(os, sp, varU, varV):

    brd=2
    dx4 = sp.dx*sp.dx*sp.dx*sp.dx

    diffU = np.zeros(varU.shape)
    diffV = np.zeros(varV.shape)

    # Diffusion along x direction:
    for i in range(brd, os.imax-1-brd):
        for j in range(brd, os.jmax-brd):
            if not os.maskU[i,j,0]:
                continue
            val = varU[i,j]
            s1 = getUV2(varU,i+1,j,val) + getUV2(varU,i,j+1,val) \
                + getUV2(varU,i-1,j,val) + getUV2(varU,i,j-1,val)
            s2 = getUV2(varU,i+1,j+1,val) + getUV2(varU,i-1,j+1,val) \
                + getUV2(varU,i-1,j-1,val) + getUV2(varU,i+1,j-1,val)
            s3 = os.maskU[i+1,j,0]*getUV2(varU,i+2,j,val) \
                    + os.maskU[i,j+1,0]*getUV2(varU,i,j+2,val) \
                    + os.maskU[i-1,j,0]*getUV2(varU,i-2,j,val) \
                    + os.maskU[i,j-1,0]*getUV2(varU,i,j-2,val)

            diffU[i,j] = sp.KBi*(20*val - 8*s1 + 2*s2 + s3)/dx4

    # Diffusion along y direction:
    for i in range(brd, os.imax-brd):
        for j in range(brd, os.jmax-1-brd):
            if not os.maskV[i,j,0]:
                continue
            val = varV[i,j]
            s1 = getUV2(varV,i+1,j,val) + getUV2(varV,i,j+1,val) \
                + getUV2(varV,i-1,j,val) + getUV2(varV,i,j-1,val)
            s2 = getUV2(varV,i+1,j+1,val) + getUV2(varV,i-1,j+1,val) \
                + getUV2(varV,i-1,j-1,val) + getUV2(varV,i+1,j-1,val)
            s3 = os.maskV[i+1,j,0]*getUV2(varV,i+2,j,val) \
                    + os.maskV[i,j+1,0]*getUV2(varV,i,j+2,val) \
                    + os.maskV[i-1,j,0]*getUV2(varV,i-2,j,val) \
                    + os.maskV[i,j-1,0]*getUV2(varV,i,j-2,val)

            diffV[i,j] = sp.KBi*(20*val - 8*s1 + 2*s2 + s3)/dx4

    return (diffU, diffV)