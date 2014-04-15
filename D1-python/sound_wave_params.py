import numpy as np
import math
from constants import gamma


#functiontype = 'sine'
functiontype = 'defined'

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

#init_functions_generation = [{'csSign':1, 'A': A}, {'csSign':-1, 'A': 0.5*A}] #SUPERPOSITION wave travelling right with amp A and left with amp 0.5 * A
#init_functions_generation = [{'csSign':-1, 'A': A}, {'csSign':1, 'A': 0.5*A}] #SUPERPOSITION wave travelling left with amp A and right with amp 0.5 * A
#init_functions_generation = [{'csSign':-1, 'A': A}, {'csSign':1, 'A': A}] #SUPERPOSITION wave travelling left with amp A and right with amp  A
#init_functions_generation = [{'csSign':-1, 'A': A}] #travelling left
init_functions_generation = [{'csSign':1, 'A': A}] #travelling right



#The following are used for the notifier(taken from notifier_params as they were wavesound specific)
plotPresCurve = False
plotRhoCurve = False
plotVelCurve = False
markPoints = True
















