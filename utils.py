import math
import numpy as np


#Get the indexed value from the matrix M. If the index is outside the
# edge, or the value is NaN, return the default value.
def getUV(M, i, j, k, defVal):
   sm = M.shape
   if i<0 or i>=sm[0] or j<0 or j>=sm[1] or k<0 or k>=sm[2]:
       return defVal
   else:
       val = M[i,j,k]
       if math.isnan(val):
           return defVal
       else:
           return val

#Get the indexed value from the 2D matrix M. If the index is outside the
# edge, or the value is NaN, return the default value.
def getUV2(M, i, j, defVal):
   sm = M.shape
   if i<0 or i>=sm[0] or j<0 or j>=sm[1]:
       return defVal
   else:
       val = M[i,j]
       if math.isnan(val):
           return defVal
       else:
           return val

# Calculate the mean value of those elements of the tuple vals that are not nan.
# If all are nan, return the default value.
def calcMean(vals, defaultVal=0):
    count = 0
    sum = 0
    for v in vals:
        if not math.isnan(v):
            count = count+1
            sum = sum+v
    if count > 0:
        return sum/count
    else:
        return defaultVal

# Smooth a matrix by a diffusion process:
def smooth(m, D=0.1, repeats=1):
    a = np.empty(m.shape)
    a[...] = m[...]
    for rep in range(0,repeats):
        for i in range(1,m.shape[0]-1):
            for j in range(1,m.shape[1]-1):
                SS = [0, 0, 0, 0]
                if ~math.isnan(m[i,j]):
                    if ~math.isnan(m[i-1,j]):
                        SS[0] = D*(m[i-1,j]-m[i,j])
                    if ~math.isnan(m[i+1,j]):
                        SS[1] = D*(m[i+1,j]-m[i,j])
                    if ~math.isnan(m[i,j-1]):
                        SS[2] = D*(m[i,j-1]-m[i,j])
                    if ~math.isnan(m[i,j+1]):
                        SS[3] = D*(m[i,j+1]-m[i,j])
                    a[i,j] = a[i,j] + sum(SS)
    return a


def maxmod(a,b):
    if a*b < 0:
        return 0
    if abs(a) > abs(b):
        return a
    else:
        return b

def minmod(a,b):
    if a*b < 0:
        return 0
    if abs(a) < abs(b):
        return a
    else:
        return b



def vertAverageUBVB(os,i,j):
    sumU = 0
    sumV = 0
    lU = 0
    lV = 0
    for k in range(0, os.kmax):
        if os.maskU[i,j,k]:
            sumU = sumU + os.UB[i,j,k]*os.DWD[i,j,k]
            lU = lU+os.DWD[i,j,k]
        if os.maskV[i,j,k]:
            sumV = sumV + os.VB[i,j,k]*os.DSD[i,j,k]
            lV = lV+os.DSD[i,j,k]
    if lU>0:
        sumU = sumU/lU
    if lV>0:
        sumV = sumV/lV
    return sumU, sumV


# Check the edges of a model domain, and make sure that no
# wet points at the border have a dry point one or two cells inside of it.
# Where this occurs, the depth at the edge is reduced to match
# the cell inside.
# The reason for this is to avoid errors when doing certain operations
# along the boundary, such as copying the values inside the boundary out
# to the boundary to get a floating boundary.
def adaptDepthField(os):
    anyFound = False
    # Left and right edges:
    for j in range(1,os.jmax-1):
        # Left:
        if os.kmm[0,j] > min(os.kmm[1,j], os.kmm[2,j]):
            os.depth[0,j] = min(os.depth[1,j], os.depth[2,j])
            os.depth[1,j] = os.depth[0,j]
            anyFound = True
        # Right:
        if os.kmm[-1,j] > min(os.kmm[-2,j],os.kmm[-3,j]):
            os.depth[-1,j] = min(os.depth[-2,j], os.depth[-3,j])
            os.depth[-2,j] = os.depth[-1,j]
            anyFound = True
    for i in range(1, os.imax-1):
        # Bottom:
        if os.kmm[i,0] > min(os.kmm[i,1],os.kmm[i,2]):
            os.depth[i,0] = min(os.depth[i,1], os.depth[i,2])
            os.depth[i,1] = os.depth[i,0]
            anyFound = True
        # Top:
        if os.kmm[i,-1] > min(os.kmm[i,-2],os.kmm[i,-3]):
            os.depth[i,-1] = min(os.depth[i,-2],os.depth[i,-3])
            os.depth[i,-2] = os.depth[i,-1]
            anyFound = True

    if anyFound:
        os.calcKmmDzz()


