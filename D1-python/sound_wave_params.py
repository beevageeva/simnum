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

#init_functions_generation = [{'csSign':1, 'A': A}, {'csSign':-1, 'A': 0.5*A}] #SUPERPOSITION wave travelling right with amp A and left with amp 0.5 * A
#init_functions_generation = [{'csSign':-1, 'A': A}, {'csSign':1, 'A': 0.5*A}] #SUPERPOSITION wave travelling left with amp A and right with amp 0.5 * A
#init_functions_generation = [{'csSign':-1, 'A': A}, {'csSign':1, 'A': A}] #SUPERPOSITION wave travelling left with amp A and right with amp  A
#init_functions_generation = [{'csSign':-1, 'A': A}] #travelling left
init_functions_generation = [{'csSign':1, 'A': A}] #travelling right
