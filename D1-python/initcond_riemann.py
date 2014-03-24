import numpy as np
import sys,math
from constants import gamma
import riemann_params



def getInitialPresRhoVel(z):
	#check (rhoLeft > rhoRight and presLeft > presRight) or (rhoLeft < rhoRight and presLeft < presRight) otherwise initial conditions are invalid
	if  (riemann_params.rhoLeft > riemann_params.rhoRight and riemann_params.presLeft < riemann_params.presRight) or  (riemann_params.rhoLeft < riemann_params.rhoRight and riemann_params.presLeft > riemann_params.presRight):
		print("invalid initial conditions : (rhoLeft >= rhoRight and presLeft >= presRight) or (rhoLeft <= rhoRight and presLeft <= presRight)")
		sys.exit(0)
	#check rhoLeft > rhoRight or inverse them
	if(riemann_params.rhoLeft < riemann_params.rhoRight):
		print("!!!!!!Rholeft>RhoRight INVERSE")
		temp = riemann_params.rhoRight
		riemann_params.rhoRight = riemann_params.rhoLeft
		riemann_params.rhoLeft = temp
		temp = riemann_params.presRight
		riemann_params.presRight = riemann_params.presLeft
		riemann_params.presLeft = temp
		temp = riemann_params.velRight
		riemann_params.velRight = riemann_params.velLeft
		riemann_params.velLeft = temp
		
	from common import getDz
	#smooth functions
	ww = 10.0 * getDz()
	pres0 =np.add(np.multiply((0.5 *  (riemann_params.presLeft - riemann_params.presRight)),np.subtract(1.0, np.tanh(np.divide((np.subtract(z,riemann_params.zC)),ww)))),riemann_params.presRight)
	rho0 =np.add(np.multiply((0.5 *  (riemann_params.rhoLeft - riemann_params.rhoRight)),np.subtract(1.0, np.tanh(np.divide((np.subtract(z,riemann_params.zC)),ww)))),riemann_params.rhoRight)
	vel0 =np.add(np.multiply((0.5 *  (riemann_params.velLeft - riemann_params.velRight)),np.subtract(1.0, np.tanh(np.divide((np.subtract(z,riemann_params.zC)),ww)))),riemann_params.velRight)
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
			

def getCsLeft():
	return math.sqrt(gamma * riemann_params.presLeft / riemann_params.rhoLeft)

def getCsRight():
	return math.sqrt(gamma * riemann_params.presRight / riemann_params.rhoRight)

