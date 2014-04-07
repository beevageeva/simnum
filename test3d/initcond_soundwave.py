import numpy as np
import sys,math
from constants import gamma
from sound_wave_params import functiontype




def getV00():
	from sound_wave_params import v00
	return v00

def getP00():
	from sound_wave_params import p00
	return p00

def getRho00():
	from sound_wave_params import rho00
	return rho00

def getCs0():
	from sound_wave_params import p00, rho00
	cs = math.sqrt(gamma *  p00 / rho00)
	return cs


def getInitialPresRhoVel(z):
	from sound_wave_params import A, p00, rho00, v00
	cs00 = math.sqrt(gamma * p00 / rho00)
	arg = z[0] #only first direction	
	if functiontype == "sine":
		from sound_wave_sine_params import wl, phi
		r =  A * np.sin(np.multiply((2.0 * math.pi/wl),arg) + phi )
	elif functiontype == "defined":
		from sound_wave_defined_params import w
		r = w(arg)
	v1 = v00 + cs00 *  r
	vel = np.dstack((v1,np.zeros(v1.shape)))
	return {'pres': p00 + gamma * p00 * r  , 'rho': rho00 + rho00 * r , 'vel': vel } 

			

