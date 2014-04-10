import numpy as np
import math
from constants import gamma


functiontype = 'sine'
#functiontype = 'defined'

rho00 = 1.0
p00 = 1.0


v00 = 0.0

#cs00 = math.sqrt(gamma * p00 / rho00)
#v00 = - cs00 /  5.5
#v00 =  0.5 * cs00
#v00 = - cs00
#v00 =  cs00

A = 3.0 * 10.0 ** (-4)
#A = 5.0 * 10.0 ** (-2)
#A = 10.0 ** (-2) #this will work fine as well with fg scheme

#periodic = "y" # can be "x", "y"  or "r" (periodic), "d"
#periodic = "r" 
#periodic = "x"
periodic = "d1" 






