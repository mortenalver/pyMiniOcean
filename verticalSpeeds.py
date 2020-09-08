from utils import getUV

# Compute vertical speeds, from bottom up per horizontal location. We
# already know the horizontal current speeds, and based on the assumption
# of incompressibility, we require the flow out of a cell to equal the flow
# into it at all times. Going from bottom to top, we find the vertical
# speed at the top of the cell that balances the flows through all the
# other edges. At the bottom, there is 0 flow though the bottom of the
# cell, and for cells further up the water column we have already computed
# the flow in the previous iteration.
# The exception to the flow balance requirement is the surface layer, whose 
# cells are allowed to have variable volume - this is handled in the same way 
# as the other cells except that the vertical speed at the top represents the 
# rate of change of the elevation instead of a flow through the top edge.
def computeVerticalSpeeds(os, sp):
	for i in range(1,os.imax-1):
		for j in range(1,os.jmax-1):
			os.W[i,j,os.kmm[i,j]] = 0 # No current through bottom
			for k in range(os.kmm[i,j]-1, -1, -1):
				# Current through upper boundary is calculated to balance
				# flow rates through the other ones.

				# Calculate mean cell heights on both borders in x and y
				# direction, to find areas of cell interfaces:
				meanHeightU = (0.5*(os.cellHeights[i-1,j,k] + os.cellHeights[i,j,k]), 0.5*(os.cellHeights[i,j,k] + os.cellHeights[i+1,j,k]))
				meanHeightV = (0.5*(os.cellHeights[i,j-1,k] + os.cellHeights[i,j,k]),  0.5*(os.cellHeights[i,j,k] + os.cellHeights[i,j+1,k]))
				
				flowDiff = os.W[i,j,k+1]*sp.dx*sp.dx \
					+ getUV(os.U,i-1,j,k,0)*sp.dx*meanHeightU[0] - getUV(os.U,i,j,k,0)*sp.dx*meanHeightU[1] \
					+ getUV(os.V,i,j-1,k,0)*sp.dx*meanHeightV[0] - getUV(os.V,i,j,k,0)*sp.dx*meanHeightV[1]
				os.W[i,j,k] = flowDiff/(sp.dx*sp.dx)


				#if i==15 and j==4 and k==0:
				#	1

				
