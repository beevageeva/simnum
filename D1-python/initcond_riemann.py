import numpy as np
import sys,math
from constants import gamma



def getInitialPresRhoVel(z):
	from riemann_params import presLeft,rhoLeft,velLeft,presRight,zC,rhoRight,velRight
	from common import getDz
	#smooth functions
	ww = 10.0 * getDz()
	pres0 =np.add(np.multiply((0.5 *  (presLeft - presRight)),np.subtract(1.0, np.tanh(np.divide((np.subtract(z,zC)),ww)))), presRight)
	rho0 =np.add(np.multiply((0.5 *  (rhoLeft - rhoRight)),np.subtract(1.0, np.tanh(np.divide((np.subtract(z,zC)),ww)))), rhoRight)
	vel0 =np.add(np.multiply((0.5 *  (velLeft - velRight)),np.subtract(1.0, np.tanh(np.divide((np.subtract(z,zC)),ww)))), velRight)
#	n = z.shape[0]
#	pres0 = np.zeros(n)
#	rho0 = np.zeros(n)
#	vel0 = np.zeros(n)
#	for i in range(0,len(z)):
#		if(z[i] <=zC):
#			pres0[i] = presLeft
#			rho0[i] = rhoLeft
#			vel0[i] = velLeft
#		else:
#			pres0[i] = presRight
#			rho0[i] = rhoRight
#			vel0[i] = velRight
	return {'pres': pres0  , 'rho': rho0 , 'vel': vel0 } 
			

