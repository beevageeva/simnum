notifierType = "visual"
nstepsPlot = 1 

#projections = None
from constants import z0_0, zf_0, z0_1, zf_1

#second element in array is whether to calculate maxpoint velocity: I don't know what happens when multiple points for the max
projections = {"dim0": [0.5 * (zf_0 + z0_0), False], "dim1": [0.5 * (zf_1 + z0_1), False], "color": True}


#use it in wave packet
plotVelFFT = False
#plotVelFFT = True

#plotAnalitical = True
plotAnalitical = False
