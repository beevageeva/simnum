import numpy as np
import math
from constants import gamma


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

#argType = "y" # can be "x", "y"  or "r", "d1"


#these sholuld be defined only in case of r argType: see below
#velDir = "x"
#velDir = "y"
#velDir = "d1"
#argType = "x"
argType = "y"
#argType = "r"
#argType = "d1" 
periodicType = "repeat" # repeat| diff



#be sure 
if(argType == "d1"):
	velDir = "d1"

elif(argType == "x"):
	velDir = "x"

elif(argType == "y"):
	velDir = "y"




