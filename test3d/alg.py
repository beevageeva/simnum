from constants import gamma
import numpy as np


def power2Vec(vfr):
	return vfr[:,:,0] ** 2 + vfr[:,:,1] ** 2


def getInitialUcUe(rho, v , p):
	uc = np.dstack((rho * v[:,:,0], rho * v[:,:,1] )) #cannot directly multiply 2d and 3d array
	print("uc shapein getInitialUcUe")
	print(uc.shape)
	print("UC ini")
	print(uc)
	ue = np.add(np.divide(p,(gamma - 1.0)),0.5 * np.multiply(rho, power2Vec(v)))
	return {'uc': uc, 'ue': ue}

def recalculateVelPres(rho, uc, ue):
	print("Rho=")
	print(rho)
	print("Uc in recalculate Vel Pres ")
	print(uc)

	v = np.dstack((np.divide(uc[:,:,0], rho ), np.divide(uc[:,:,1], rho )  )) #cannot directly divide 3d to 2d array
	t1 = np.subtract(ue,np.divide(power2Vec(uc), np.multiply(rho, 2.0)))	
	p = np.multiply((gamma - 1.0), t1)
	print("vel ")
	print(v)
	return {'vel': v, 'pres': p}

def recalculateFluxes(rho, uc, ue, v, p):
	fm = uc
	fe1 = ue + p	
	fe = np.dstack((fe1 * v[:,:,0], fe1 * v[:,:,1] )) #cannot directly multiply 2d and 3d array
	fc1 = np.multiply(rho, v[:,:,0] ** 2) + p
	fc2 = np.multiply(rho, np.multiply(v[:,:,0], v[:,:,1])) #Fc_xz = Fc_zx store only ince the value
	fc3 = np.multiply(rho, v[:,:,1] ** 2) + p
	fc = np.dstack((fc1,fc2,fc3))
	return {'fm': fm, 'fc': fc, 'fe': fe}

def getTimestep(v, p, rho):
	from constants import fcfl, verbose
	from common import getDz
	dz = getDz()
	t1 =  np.divide(p, rho)
	if(np.any(t1<0)):
		print("tehere is p/rho < 0 in getTimestep")
		print("p<0")
		print(p[p<0])
		print("rho<0")
		print(rho[rho<0])

		return 0	
	cs = np.sqrt(gamma *  t1)
	velMod = np.sqrt(power2Vec(v))	
	smax = np.max(np.concatenate([np.absolute(velMod + cs), np.absolute(velMod - cs)]))
	dt = float(dz * fcfl ) / smax
	return dt


def lrBoundaryConditions(array, skip=0):
	from constants import problemType
	n = array.shape[0] - 1
	if problemType == "riemann":
		array = np.insert(array, 0, array[0,:], axis = 0)
		array = np.insert(array, n, array[n,:], axis = 0)
		array = np.insert(array, 0, array[:,0], axis = 1)
		array = np.insert(array, n, array[:,n], axis = 1)
		return array
	elif problemType == "soundwave":
		array = np.insert(array, 0, array[n-skip,:], axis = 0)
		array = np.insert(array, n, array[1+skip,:], axis = 0)
		array = np.insert(array, 0, array[:,n-skip], axis = 1)
		array = np.insert(array, n, array[:,1+skip], axis = 1)
		return array
	else:
		print("problemtype %s  not implemented " % problemType)


from constants import schemeType

if schemeType == "fg":
	
	def calcIntermU(u, f, dt):
		print("calcIntermU")
		from common import getDz
		from constants import nint
		lambdaParam = dt / getDz()
		if(u.ndim == 2): #case of um and ue that are not vectors
			res = np.zeros((nint+1, nint+1))
			#res = np.array(nint+1, nint+1) #TODO do not initialize how?
		else:
			res = np.zeros((nint+1, nint+1, 2))
			#res = np.array(nint+1, nint+1, 2)
		for i in range(1, nint+2):
			for j in range(1, nint+2):
				#points displaced right +1 
				if(u.ndim == 2): #case of um and ue that are not vectors
					val = 0.25 * (u[i][j] + u[i-1][j-1] + u[i][j-1] + u[i-1][j]) - 0.25 * lambdaParam  * ((f[i][j][0] - f[i-1][j][0] + f[i][j-1][0] - f[i-1][j-1][0]) + (f[i][j][1] - f[i][j-1][1] + f[i-1][j][1] - f[i-1][j-1][1] ))
					res[i-1][j-1] = val
				else:
					val = 0.25 * (u[i][j][0] + u[i-1][j-1][0] + u[i-1][j][0] + u[i][j-1][0]) - 0.25 * lambdaParam  * ((f[i][j][0] - f[i-1][j][0] + f[i][j-1][0] - f[i-1][j-1][0]) + (f[i][j][1] - f[i][j-1][1] + f[i-1][j][1] - f[i-1][j-1][1]))
					val2 = 0.25 * (u[i][j][1] + u[i-1][j-1][1] + u[i][j-1][1] + u[i-1][j][1]) - 0.25 * lambdaParam  * ((f[i][j][1] - f[i-1][j][1] + f[i][j-1][1] - f[i-1][j-1][1]) + (f[i][j][2] - f[i][j-1][2] + f[i-1][j][2] - f[i-1][j-1][2] ))
					res[i-1][j-1][0] = val
					res[i-1][j-1][1] = val2
					#res[i-1].append([val,val2])
	
		#left and right boundary condition  skip one point !!! both right and left the intermediate array will have nint + 3 points see array limits
		res = lrBoundaryConditions(res, 1)
		return res
	
	
	def calcFinalU(u, intermF, dt):
		print("calcFinalU")
		from common import getDz
		from constants import nint
		lambdaParam = dt / getDz()
		if(u.ndim == 2): #case of um and ue that are not vectors
			res = np.zeros((nint+2, nint+2))
			#res = np.array(nint+1, nint+1) #TODO do not initialize how?
		else:
			res = np.zeros((nint+2, nint+2, 2))
			#res = np.array(nint+1, nint+1, 2)
		for i in range(0, nint+2):
			for j in range(0, nint+2):
				if(u.ndim == 2): #case of um and ue that are not vectors
					val = u[i][j] - lambdaParam  * 0.5 * ((intermF[i+1][j][0] - intermF[i][j][0] + intermF[i+1][j+1][0] - intermF[i][j+1][0]) + (intermF[i][j+1][1] - intermF[i][j][1] + intermF[i+1][j+1][1] - intermF[i+1][j][1]))
					res[i][j]=val
				else:
					val = u[i][j][0] - lambdaParam  * 0.5 * ((intermF[i+1][j][0] - intermF[i][j][0] + intermF[i+1][j+1][0] - intermF[i][j+1][0]) + (intermF[i][j+1][1] - intermF[i][j][1] + intermF[i+1][j+1][1] - intermF[i+1][j][1]))
					val2 = u[i][j][1] - lambdaParam  * 0.5 * ((intermF[i+1][j][1] - intermF[i][j][1]+intermF[i+1][j+1][1] - intermF[i][j+1][1]) + (intermF[i][j+1][2] - intermF[i][j][2]+ intermF[i+1][j+1][2] - intermF[i+1][j][2]))
					res[i][j][0] = val
					res[i][j][1] = val2
		#no more boundary conditions because intermediate array alreday has nint + 3 points
		return np.array(res)
		
	
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
		if(u.ndim == 2): #case of um and ue that are not vectors
			res = np.zeros((nint, nint))
		else:
			res = np.zeros((nint, nint, 2))
		for i in range(1, nint+1):
			for j in range(1, nint+1):
				if(u.ndim == 2): #case of um and ue that are not vectors
					val = 0.25 * (u[i-1][j-1] + u[i-1][j+1] + u[i+1][j+1] + u[i+1][j-1]) - 0.25 * lambdaParam  * ((f[i+1][j-1][0] - f[i-1][j-1][0]+ f[i+1][j+1][0] - f[i-1][j+1][0]) + (f[i-1][j+1][1] - f[i-1][j-1][1] + f[i+1][j+1][1] - f[i+1][j-1][1])) 
					res[i-1][j-1] = val
				else:
					val = 0.25 * (u[i-1][j-1][0] + u[i-1][j+1][0] + u[i+1][j+1][0] + u[i+1][j-1][0]) - 0.25 * lambdaParam  * ((f[i+1][j-1][0] - f[i-1][j-1][0] + f[i+1][j+1][0] - f[i-1][j+1][0]) + (f[i-1][j+1][1] - f[i-1][j-1][1] + f[i+1][j+1][1] - f[i+1][j-1][1])) 
					val2 = 0.25 * (u[i-1][j-1][1] + u[i-1][j+1][1] + u[i+1][j+1][1] + u[i+1][j-1][1]) - 0.25 * lambdaParam  * ((f[i+1][j-1][1] - f[i-1][j-1][1]+ f[i+1][j+1][1] - f[i-1][j+1][1]) + (f[i+1][j+1][2] - f[i+1][j-1][2]+ f[i-1][j+1][2] - f[i-1][j-1][2])) 
					res[i-1][j-1][0] = val
					res[i-1][j-1][1] = val2
		print("res before BC")
		print(res)
		res = lrBoundaryConditions(res)
		print("res AFTER BC")
		print(res)
		return res
	
	def recalculateU(rho, uc, ue, fm , fc, fe, dt):
		finalRho = calcSingleStepU(rho, fm, dt)
		print("uc  (before recalculateU) = ")
		print(uc)	
		print("fc (recalculateU) = ")
		print(fc)	
		finalUc = calcSingleStepU(uc, fc, dt)	
		print("uc  (after recalculateU) = ")
		print(finalUc)	
		finalUe = calcSingleStepU(ue, fe, dt)
		return {"rho": finalRho, "uc": finalUc, "ue": finalUe}

else:
	puts("Scheme type not implemented %s" % schemeType)



