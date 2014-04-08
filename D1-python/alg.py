from constants import gamma
import numpy as np

def getInitialUcUe(rho, v , p):
	uc = np.multiply(rho, v)
	ue = np.add(np.divide(p,(gamma - 1.0)),0.5 * np.multiply(rho, np.power(v, 2.0)))
	print("getInitialUcUe uc:")
	print(" ".join(map(str, uc)))
	print("getInitialUcUe ue:")
	print(" ".join(map(str, ue)))
	
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
	print("getTimestep pres")
	print(" ".join(map(str, p)))
	print("getTimestep rho")
	print(" ".join(map(str, rho)))
	print("getTimestep vel")
	print(" ".join(map(str, v)))
	from constants import fcfl, verbose
	from common import getDz
	dz = getDz()
	t1 =  np.divide(p, rho)
	if(np.any(t1<0)):
		return 0	
	cs = np.sqrt(gamma *  t1)
	print("getTimestep cs")
	print(" ".join(map(str, cs)))
	smax = np.max(np.concatenate([np.absolute(v + cs), np.absolute(v - cs)]))
	dt = float(dz * fcfl ) / smax
	print("getTimestep %E" % dt)	
	return dt


def lrBoundaryConditions(array, skip=0):
	from constants import problemType
	if problemType == "riemann":
		from initcond_riemann import lrBoundaryConditions as bc
	elif problemType == "soundwave":
		from initcond_soundwave import lrBoundaryConditions as bc
	else:
		print("problemtype %s  not implemented " % problemType)
		return
	bc(array, skip)

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
			print(val)
			res.append(val)
		#left and right boundary condition  skip one point !!! both right and left the intermediate array will have nint + 3 points see array limits
		print("calcIntermStep before bc")
		print(res)	
		lrBoundaryConditions(res, 1)
		print("calcIntermStep after bc")
		print(res)	
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
		print("calcIntermRho ")
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
			print(val)	
			res.append(val)
		print("calcSingleStep before bc")
		print(res)	
		lrBoundaryConditions(res)
		print("calcSingleStep after bc")
		print(res)	
		return res
	
	def recalculateU(rho, uc, ue, fm , fc, fe, dt):
		print("recalculateRho")
		finalRho = calcSingleStepU(rho, fm, dt)	
		print("recalculateUc")
		finalUc = calcSingleStepU(uc, fc, dt)	
		print("recalculateUe")
		finalUe = calcSingleStepU(ue, fe, dt)
		return {"rho": finalRho, "uc": finalUc, "ue": finalUe}

else:
	puts("Scheme type not implemented %s" % schemeType)



