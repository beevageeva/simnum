import numpy as np
import sys
from constants import gamma
from base_model import BaseModel
from math import sqrt

from riemann_params import riemann_problemType


mark4Points = True


def getCs(p, rho):
	t1 = p / rho
	if t1 < 0:
		return 0	
	return sqrt(gamma * t1)	

def getCsRWave():
	from initcond_riemann import getCsLeft
	from riemann_params import velLeft
	return (velLeft-getCsLeft()) 


def checkExpr14(cs1, cs2, pJ):
	print("-----------check expr 5.14----------")
	if(cs1==0):
		print("cs1=0")
		return
	t1 = (cs2 * cs2)/ (cs1 * cs1)
	gd = (gamma + 1.0) / (gamma - 1.0)
	t2 = pJ * (gd + pJ) / (1.0 + gd * pJ)
	print("t1=%E, t2=%E, absdif=%E" % (t1, t2, abs(t1-t2)))

def checkExpr15(v2, v1, cs1, pJ):
	print("-----------check expr 5.15----------")
	t23 = (gamma + 1.0)*(pJ -1.0)/ (2.0 * gamma) + 1.0
	if(t23<0):
		print("negative for sqrt")
		return
	t22 = sqrt(t23)
	t21 = (pJ - 1.0) / t22
	t2 = v1 + (cs1 / gamma) * t21
	print("t1=%E, t2=%E, absdif=%E" % (v2, t2, abs(v2-t2)))
	 

def checkExpr16(v3, v4, cs4, p3, p4):
	print("-----------check expr 5.16----------")
	t21 = 1.0 - (p3/p4) ** ((gamma-1.0)/(2.0*gamma)) 
	t2 = v4 + (2.0 * cs4 * t21)/(gamma - 1.0)
	print("t1=%E, t2=%E, absdif=%E" % (v3, t2, abs(v3-t2)))
	
def checkExpr17(p4, p1, pJ, cs4, v4, v1, cs1):
	print("-----------check expr 5.17----------")
	t1 = p4 / p1
	t25 = ((gamma + 1.0) / (2.0 * gamma)) * (pJ - 1.0) + 1.0	
	if(t25<0):
		print("negative for sqrt")
		return
	t24 = sqrt(t25)
	t23 = (pJ - 1.0) / t24
	t22 = v4 - v1 - (cs1/gamma) * t23
	t21 = 1.0 + ((gamma - 1.0)/ (2.0 * cs4)) * t22 
	t2 = pJ * t21 ** ((-2.0 * gamma)/ (gamma -1.0))
	print("t1=%E, t2=%E, absdif=%E" % (t1, t2, abs(t1-t2)))




class Model(BaseModel):
	
	def __init__(self):
		BaseModel.__init__(self)



	def checkExpressions(self):
		#calculate delta points before and after the points zC , rwPoint, shPoint
		n = len(self.pres)
		delta = int(n / 40)
		from riemann_params import  presLeft, presRight, rhoLeft, rhoRight, velLeft, velRight
		from common import getZIndex
		zIndexSh = getZIndex(self.shPoint) 		
		zIndexZc = getZIndex(self.zcPoint) 		
		zIndexRw = getZIndex(self.rwPoint) 	
		p1Index = zIndexSh + delta
		if(p1Index >=n):
			p1Index = n-1
		p2Index = zIndexSh - delta
		p3Index = zIndexZc + delta
		p4Index = zIndexRw - delta
		if(p4Index < 0):
			p4Index = 0
		
		#keep the points 
		if(mark4Points):
			self.point1 = self.z[p1Index]		
			self.point2 = self.z[p2Index]		
			self.point3 = self.z[p3Index]		
			self.point4 = self.z[p4Index]		
		
		print("--------------delta=%d--------------------------" % delta)	
		print("pres1=%E" % self.pres[p1Index])
		print("rho1=%E" % self.rho[p1Index])
		print("vel1=%E" % self.vel[p1Index])
		cs1 = getCs(self.pres[p1Index], self.rho[p1Index])
		print("cs1=%E" % cs1)
		print("pres2=%E" % self.pres[p2Index])
		print("rho2=%E" % self.rho[p2Index])
		print("vel2=%E" % self.vel[p2Index])
		cs2 = getCs(self.pres[p2Index], self.rho[p2Index])
		print("cs2=%E" % cs2)
		#uncomment following to see that it is not practical to calculate csShock like this(because of discontinuity and oscillations)
		#print("csShock calculated in shPoint=%E" % getCs(self.pres[zIndexSh], self.rho[zIndexSh])
		#SHOCK ANALYTIC
		if riemann_problemType in ["shock_tube", "exp_vacuum"]:
			print("csShock=%E" % self.getCsShock())
		print("pres3=%E" % self.pres[p3Index])
		print("rho3=%E" % self.rho[p3Index])
		print("vel3=%E" % self.vel[p3Index])
		cs3 = getCs(self.pres[p3Index], self.rho[p3Index])
		print("cs3=%E" % cs3)
		print("pres4=%E" % self.pres[p4Index])
		print("rho4=%E" % self.rho[p4Index])
		print("vel4=%E" % self.vel[p4Index])
		cs4 = getCs(self.pres[p4Index], self.rho[p4Index])
		print("cs4=%E" % cs4)
		print("csRW=%E" % getCsRWave())
		print("-------check const----")
		print("presRight=%E, p1=%E, abs(presRight - p1) = %E" % (presRight,self.pres[p1Index], abs(presRight - self.pres[p1Index])))
		print("rhoRight=%E, p1=%E, abs(rhoRight - p1) = %E" % (rhoRight,self.rho[p1Index], abs(rhoRight - self.rho[p1Index])))
		print("velRight=%E, p1=%E, abs(velRight - p1) = %E" % (velRight,self.vel[p1Index], abs(velRight - self.vel[p1Index])))
		print("presLeft=%E, p6=%E, abs(presLeft - p6) = %E" % (presLeft, self.pres[p4Index], abs(presLeft - self.pres[p4Index])))
		print("rhoLeft=%E, p6=%E, abs(rhoLeft - p6) = %E" % (rhoLeft, self.rho[p4Index], abs(rhoLeft - self.rho[p4Index])))
		print("velLeft=%E, p6=%E, abs(velLeft - p6) = %E" % (velLeft, self.vel[p4Index], abs(velLeft - self.vel[p4Index])))
		print("p2=%E, p3=%E, abs(p2 - p3) = %E" % (self.pres[p2Index],self.pres[p3Index], abs(self.pres[p2Index] - self.pres[p3Index])))
		print("v2=%E, v3=%E, abs(v2 - v3) = %E" % (self.vel[p2Index],self.vel[p3Index], abs(self.vel[p2Index] - self.vel[p3Index])))
		print("----------------------------------------")
		pJ = self.pres[p2Index] / self.pres[p1Index]	
		checkExpr14(cs1, cs2, pJ)
		#checkExpr15(v2, v1, cs1, pJ)
		checkExpr15( self.vel[p2Index], self.vel[p1Index], cs1, pJ)
		#checkExpr16(v3, v4, cs4, p3, p4)
		checkExpr16(self.vel[p3Index], self.vel[p4Index], cs4, self.pres[p3Index], self.pres[p4Index])
		#checkExpr17(p4, p1, pJ, cs4, v4, v1, cs1)
		checkExpr17(self.pres[p4Index], self.pres[p1Index], pJ, cs4, self.vel[p4Index], self.vel[p1Index], cs1)
		print("------------------end----------------------")


	def updateValuesModel(self, dt, time):
		from common import getDz, getZIndex
		from riemann_params import zC
		delta = getDz()
		from analyze_functions import getFirstIndexDifferentLeft,getFirstIndexDifferentRight, getFirstIndexConstant
		dind = getZIndex(self.zcPoint)
		zi = getFirstIndexDifferentLeft(self.pres, delta, dind)
		self.rwPoint = self.z[zi]
		zi2 = getFirstIndexDifferentRight(self.pres, delta, dind)
		self.shPoint = self.z[zi2] 
		#ANALYTICAL DCPoint
		self.dcPoint = self.getNewDcPoint(dt)
		#ANALYTICAL RWPoint
		from common import displacedPoint
		from riemann_params import  velLeft
		from initcond_riemann import getCsLeft
		from riemann_params import timeAfterAnPoints
		if(time>=timeAfterAnPoints):
			#calculate zc after the time
			self.zcPoint = self.z[getFirstIndexConstant(self.pres, zi,delta)]
			#show analytical points
			if hasattr(self, 'rwPointAn'):
				self.rwPointAn = self.getNewRwPoint(dt)
				if riemann_problemType in ["shock_tube", "exp_vacuum"]:
					self.shPointAn = self.getNewShPoint(dt)
			else:
				self.initPointsAn()
			#calculate cs in all 6 points and check expressions
			print("Update Values in Model at Time %4.3f" % time)
			self.checkExpressions()




	def mark4PointsOnGraph(self, axtitle):
		self.notifier.markPoint(axtitle, axtitle + "Point1", self.point1)
		self.notifier.markPoint(axtitle, axtitle + "Point2", self.point2)
		self.notifier.markPoint(axtitle, axtitle + "Point3", self.point3)
		self.notifier.markPoint(axtitle, axtitle + "Point4", self.point4)



	def updateValuesNotifier(self, dt, time):
		self.notifier.updateValues("rho", self.rho)
		self.notifier.updateValues("pres", self.pres)
		self.notifier.updateValues("vel", self.vel)
		#zC is also moving
		self.notifier.markPoint("rho", "zCRho", self.zcPoint)
		self.notifier.markPoint("pres", "zCPres", self.zcPoint)
		self.notifier.markPoint("vel", "zCVel", self.zcPoint)

		self.notifier.markPoint("pres", "rwPres", self.rwPoint)
		self.notifier.markPoint("pres", "shPres", self.shPoint)
		self.notifier.markPoint("rho", "rwRho", self.rwPoint)
		self.notifier.markPoint("rho", "shRho", self.shPoint)
#		#the discontinuity point in density ..
		self.notifier.markPoint("rho", "dcRho", self.dcPoint)
		self.notifier.markPoint("vel", "rwVel", self.rwPoint)
		self.notifier.markPoint("vel", "shVel", self.shPoint)
		
		if hasattr(self, 'rwPointAn'): #I could also check by time as in updateValuesModel
			self.notifier.markPoint("pres", "rwPresAn", self.rwPointAn)
			self.notifier.markPoint("rho", "rwRhoAn", self.rwPointAn)
			self.notifier.markPoint("vel", "rwVelAn", self.rwPointAn)
			if riemann_problemType in ["shock_tube", "exp_vacuum"]:
				self.notifier.markPoint("vel", "shVelAn", self.shPointAn)
				self.notifier.markPoint("pres", "shPresAn", self.shPointAn)
				self.notifier.markPoint("rho", "shRhoAn", self.shPointAn)
			if(mark4Points):
				self.mark4PointsOnGraph("pres")
				self.mark4PointsOnGraph("rho")
				self.mark4PointsOnGraph("vel")



	def getInitialValues(self):
		return [self.pres, self.rho, self.vel]


	#initialize rwPointAn and shPointAn from points determined emp(analyze_functions) after a time when it's better represented
	def initPointsAn(self):
		#TODO
		self.rwPointAn = self.rwPoint
		if riemann_problemType in ["shock_tube", "exp_vacuum"]:
			self.shPointAn = self.shPoint


	def additionalInit(self):
		from riemann_params import zC
		from common import getDz
		#for the moment the point name is the key in the hash - I should take into account the ax title and be able to repeat the pointNames in multiple axes - I have to have different names for different axes!

		delta = getDz()
		from analyze_functions import getFirstIndexDifferentLeft,getFirstIndexDifferentRight
		zi = getFirstIndexDifferentLeft(self.pres, delta, zC)
		self.rwPoint = self.z[zi]
		zi = getFirstIndexDifferentRight(self.pres, delta, zC)
		self.shPoint = self.z[zi]
		#dcPoint analytical, but I can't determine position at time t empirically : see analyze_functions(#TODO)
		self.dcPoint = zC
		self.zcPoint = zC 

	
	def getNewDcPoint(self, dt):
		from common import displacedPoint, getZIndex
		from riemann_params import rhoRight, zC
		#zIndex = getZIndex(self.dcPoint)	
		#take 2 points right to this one because of the discontinuity
		#zIndex +=2
		zIndex = getZIndex(0.5*(zC + self.shPoint))	
		newz = displacedPoint(self.dcPoint, self.vel[zIndex], dt)
		return newz



	if riemann_problemType in ["shock_tube", "exp_vacuum"]:

		def getCsShock(self):
			from common import  getZIndex
			from riemann_params import rhoRight
			#to be sure I have right point 	
			zIndex = getZIndex(0.5*(self.shPoint + self.dcPoint))	
			#zIndex = getZIndex(self.point2)  #point2 is created after
			v = (self.vel[zIndex] * self.rho[zIndex]) / (self.rho[zIndex] - rhoRight)
			return v

		def getNewShPoint(self, dt):
			v = self.getCsShock()
			#do not use displacedPoint from common.py as it make periodic assumption when it goes away from domain, but in this case it should remain at the end
			newz = self.shPointAn + v * dt
			from constants import zf
			if newz > zf:
				return zf
			return newz




	
	def getNewRwPoint(self, dt):
		#do not use displacedPoint from common.py as it make periodic assumption when it goes away from domain, but in this case it should remain at the end
		newz = self.rwPointAn + getCsRWave()* dt
		from constants import z0
		if newz < z0:
			return z0
		return newz


