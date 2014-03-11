import numpy as np
import math
from scipy.special import jv

rho00 = 1.0
p00 = 1.0
v00 = 0.0
A = 3.0 * 10.0 ** (-4)

def w(z):
	return jv(0,  np.multiply(z, 2.2)) * (np.add(np.tanh((np.subtract(z,4.3))/0.3),1.0)) * (np.subtract(1.0,np.tanh((np.subtract(z,6.2))/0.3))) 
