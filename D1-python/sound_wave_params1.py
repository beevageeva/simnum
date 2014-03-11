import numpy as np
import math
from constants import gamma

rho00 = 1.0
p00 = 1.0
v00 = 0.0
#v00 = - (gamma * p00) / (rho00 * 5.5)
#v00 =  (0.5 * gamma * p00) / rho00 
#v00 = - (0.9 * gamma * p00) / rho00 
#v00 =  - (gamma * p00) / rho00 
#v00 =  - (0.79 * gamma * p00) / rho00 
A = 3.0 * 10.0 ** (-4)
def w(z):
	from constants import z0, zf
	return np.sin( ((2.0 * math.pi)  / (zf - z0))* np.subtract(z, z0) + math.pi / 6.0)
