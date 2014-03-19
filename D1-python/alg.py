from constants import gamma
import numpy as np

def getInitialUcUe(rho, v , p):
	uc = np.multiply(rho, v)
	ue = np.add(np.divide(p,(gamma - 1.0)),0.5 * np.multiply(rho, np.power(v, 2.0)))
	return {'uc': uc, 'ue': ue}

def recalculateVelPres(rho, uc, ue):
	v = np.divide(uc, rho)
	t1 = np.subtract(ue,np.divide(np.power(uc, 2.0), np.multiply(rho, 2.0)))	
	p = np.multiply((gamma - 1.0), t1)
	return {'vel': v, 'pres': p}

def recalculateFluxes(rho, uc, ue, v, p):
	fm = uc
	fc = np.divide(np.power(uc, 2.0), rho) + p
	fe = np.multiply((ue + p), v)
	return {'fm': fm, 'fc': fc, 'fe': fe}

def getTimestep(v, p, rho):
	from constants import fcfl, verbose
	from common import getDz
	dz = getDz()
	t1 =  np.divide(p, rho)
	print("LAST 1025, pres, rho")
	print(p[1025])
	print(rho[1025])
	if(np.any(t1<0)):
		print("!!!!!!!!!!!!!!!!p/rho has elements <0 -- taking 0 -- no imaginary part(to take abs value?)")
		print(t1)
		indices = np.argwhere(t1<0).flatten()
		print("indices")
		print(indices)
		print("pres indices")
		print(p[indices])
		print("rho indices")
		print(rho[indices])
		t1[t1<0] = 0
		print("Modified t1")
		print(t1)
		import time
		time.sleep(5)	
	cs = np.sqrt(gamma *  t1)
	smax = np.max(np.concatenate([np.absolute(v + cs), np.absolute(v - cs)]))
	dt = float(dz * fcfl ) / smax
	return dt


def lrBoundaryConditions(array, skip=0):
	from constants import problemType
	n = len(array) - 1
	if problemType == "riemann":
		array.insert(0, array[0])	
		array.append(array[n])	
	elif problemType == "soundwave":
		array.insert(0, array[n - skip])
		array.append(array[1 + skip])
	else:
		print("problemtype %s  not implemented " % problemType)


from constants import schemeType

if schemeType == "fg":
	
	def calcIntermU(u, f, dt):
		from common import getDz
		from constants import nint
		lambdaParam = dt / getDz()
		res = []
		for i in range(1, nint+2):
			#points displaced right +1 
			val = 0.5 * (u[i] + u[i-1]) - 0.5 * lambdaParam  * (f[i] - f[i-1]) 
			res.append(val)
		#left and right boundary condition  skip one point !!! both right and left the intermediate array will have nint + 3 points see array limits
		lrBoundaryConditions(res, 1)
		#res.insert(0, res[nint-1])
		#res.append(res[2])
		return res
	
	
	def calcFinalU(u, intermF, dt):
		from common import getDz
		from constants import nint
		lambdaParam = dt / getDz()
		res = []
		for i in range(0, nint+2):
			val = u[i] - lambdaParam  * (intermF[i+1] - intermF[i]) 
			res.append(val)
		#no more boundary conditions because intermediate array alreday has nint + 3 points
		return res
		
	
	def recalculateU(rho, uc, ue, fm, fc ,fe, dt):
		intermRho = calcIntermU(rho, fm , dt)	
		intermUc = calcIntermU(uc, fc , dt)	
		intermUe = calcIntermU(ue, fe , dt)
		intermVelPres = recalculateVelPres(intermRho, intermUc, intermUe)
		intermVel = intermVelPres["vel"]
		intermPres = intermVelPres["pres"]
		intermFluxes = recalculateFluxes(intermRho, intermUc, intermUe, intermVel, intermPres)
		finalRho = calcFinalU(rho, intermFluxes["fm"], dt)	
		finalUc = calcFinalU(uc, intermFluxes["fc"], dt)	
		finalUe = calcFinalU(ue, intermFluxes["fe"], dt)
		return {"rho": finalRho, "uc": finalUc, "ue": finalUe}
	

elif schemeType == "lf":

	def calcSingleStepU(u,f, dt):
		from common import getDz
		from constants import nint
		lambdaParam = dt / getDz()
		res = []
		for i in range(1, nint+1):
			val = 0.5 * (u[i-1] + u[i+1]) - 0.5 * lambdaParam  * (f[i+1] - f[i-1]) 
			res.append(val)
		#periodic boundary condition: I don't see the point of defining another function
		lrBoundaryConditions(res)
		#res.insert(0, res[nint-1])
		#res.append(res[1])
		return res
	
	def recalculateU(rho, uc, ue, fm , fc, fe, dt):
		finalRho = calcSingleStepU(rho, fm, dt)	
		finalUc = calcSingleStepU(uc, fc, dt)	
		finalUe = calcSingleStepU(ue, fe, dt)
		return {"rho": finalRho, "uc": finalUc, "ue": finalUe}

else:
	puts("Scheme type not implemented %s" % schemeType)



