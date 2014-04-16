notifierType = "visual"
nstepsPlot = 10 
#I might want this set to False if I have multiple figures like in case of the sound waves when I might want to plot also the curves
#fullscreenMainfigure = True 
fullscreenMainfigure = False


#projections = None
from constants import z0_0, zf_0, z0_1, zf_1, z0_2, zf_2

#second element in array is whether to calculate maxpoint velocity: I don't know what happens when multiple points for the max
projections = {"dim0": 0.5 * (zf_0 + z0_0), "dim1": 0.5 * (zf_1 + z0_1), "dim2": 0.5 * (zf_2 + z0_2)}

