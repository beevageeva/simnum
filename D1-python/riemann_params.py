#riemann_problemType = "complete"  #this can be shock_tube, complete, exp_vacuum
riemann_problemType = "conv"  #this can be shock_tube, complete, exp_vacuum

timeAfterAnPoints = 2.0   #default
if riemann_problemType == "shock_tube":
	#shock tube
	timeAfterAnPoints = 1.0   

	#presLeft = 1.0
	#presRight = 0.1
	#rhoLeft = 1.0
	#rhoRight = 0.125
	#velLeft = 0.0
	#velRight = 0.0
	#zC = 1.5
	#reversed I will check for rhoLeft>rhoRight in initcond_riemann as I only want the case rhoLeft > rhoRight 
	presRight = 1.0
	presLeft = 0.1
	rhoRight = 1.0
	rhoLeft = 0.125
	velRight = 0.0
	velLeft = 0.0
	zC = 1.5
	
	#TEST only reverse pres INVALID pres1<pres2 and rho1>rho2
	#presLeft = 0.1
	#rhoLeft = 1.0
	#presRight = 1.0
	#rhoRight = 0.125
	#velLeft = 0.0
	#velRight = 0.0
	#zC = 1.5

elif riemann_problemType == "complete":
	#complete riemann problem example : velocities !=0
	timeAfterAnPoints = 0.004   #see readme

	presLeft = 10**5
	presRight = 10**4
	rhoLeft = 1.0
	rhoRight = 0.125
	velLeft = 100.0
	velRight = -50.0
	##velRight = -300 #greater than csShock
	##velRight = 50.0 #same sign
	##velRight = 1000.0 #same sign greater than csShock

	#reverse signs strange things 
	#velLeft = -100
	##velRight = 50
	#velRight = 10000 #

	zC = 0.0


elif riemann_problemType == "exp_vacuum":
	timeAfterAnPoints = 0.6   #used for ..see readme
	#expansion vacuum HT
	presLeft = 1.0
	rhoLeft = 1.0
	rhoRight = rhoLeft * 7 * 10 ** (-3)
	#case a	
	#presRight = 0.1 
	#case b
	presRight = rhoRight * 0.8
	velLeft = 0.0
	velRight = 0.0
	zC = 1.5

elif riemann_problemType == "conv":
	#converging flow
	presLeft = 1.0
	presRight = 1.0
	rhoLeft = 1.0
	rhoRight = 1.0
	velLeft = 10.0
	velRight = -10.0
	zC = 2.5
