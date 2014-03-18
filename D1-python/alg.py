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
	cs = np.sqrt(gamma *  np.divide(p, rho))
	smax = np.max(np.concatenate([np.absolute(v + cs), np.absolute(v - cs)]))
	dt = float(dz * fcfl ) / smax
	return dt

from constants import schemeType

if schemeType == "fg":
	
	def calcIntermU(u, f, dt):
		from common import getDz
		from constants import nint
		lambdaParam = dt / getDz()
		res = []
		for i in range(0, nint+1):
			val = 0.5 * (u[i] + u[i+1]) - 0.5 * lambdaParam  * (f[i+1] - f[i]) 
			res.append(val)
		#right boundary condition
		res.append(res[0])
		return res
	
	
	def calcFinalU(u, intermF, dt):
		from common import getDz
		from constants import nint
		lambdaParam = dt / getDz()
		res = []
		for i in range(0, nint+1):
			val = u[i] - lambdaParam  * (intermF[i] - intermF[i+1]) 
			res.append(val)
		#left boundary condition
		res.insert(0, res[nint])
		return res
		
	
	def recalculateU(rho, uc, ue, fm, fc ,fe, dt):
		intermRho = calcIntermU(rho, fm , dt)	
		intermUc = calcIntermU(uc, fc , dt)	
		intermUe = calcIntermU(ue, fe , dt)	
		intermVelPres = recalculateVelPres(intermRho, intermUc, intermUe)
		intermVel = intermVelPres["vel"]
		intermPres = intermVelPres["pres"]
		intermFluxes = recalculateFluxes(intermRho, intermUc, intermUe, intermVel, intermPres)
		finalRho = calcFinalU(rho, intermFluxes["fm"])	
		finalUc = calcFinalU(uc, intermFluxes["fc"])	
		finalUe = calcFinalU(ue, intermFluxes["fe"])	
		return {"rho": finalRho, "uc": finalUc, "ue": finalUe}
	

elif schemeType == "lf":

	def calcSingleStepU(u,f, dt):
		from common import getDz
		from constants import nint
		lambdaParam = dt / getDz()
		res = []
		for i in range(1, nint+2):
			val = 0.5 * (u[i-1] + u[i+1]) - 0.5 * lambdaParam  * (f[i+1] - f[i-1]) 
			res.append(val)
		#periodic boundary condition: I don't see the point of defining another function
		res.insert(0, res[nint-1])
		res.append(res[1])
		return res
	
	def recalculateU(rho, uc, ue, fm , fc, fe, dt):
		finalRho = calcSingleStepU(rho, fm, dt)	
		finalUc = calcSingleStepU(uc, fc, dt)	
		finalUe = calcSingleStepU(ue, fe, dt)
		return {"rho": finalRho, "uc": finalUc, "ue": finalUe}

else:
	puts("Scheme type not implemented %s" % schemeType)



