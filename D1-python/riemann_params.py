riemann_problemType = "shock_tube"  #this can be shock_tube, complete, exp_vacuum

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
	#u1<0(one), u2>0(one)  #see u1 and u2 def in initcond_riemann.py
	velLeft = 100.0
	velRight = -50.0
	#velRight = -300 #greater than csShock ??, but there is one u2>0 (shock wave exists)
	zC = 0.0

	#u1<0(both), u2>0(one)
	#velLeft = -1000.0
	#velRight = -50.0
	#zC = 2.0
	

	#u1<0(both), u2>0(both)
	#velLeft = -1000.0
	#velRight = 1000.0

	#the following cases are invalid because there is no rarefaction wave or no shock wave
	#u1<0 (none), u2 > 0(one) no rarefaction wave
	#velLeft = 1000.0
	#velRight = -50.0
	
	#u1<0 (one), u2 > 0(none) #no shock wave
	#velLeft = 100.0
	#velRight = -1000.0

	#u1<0 (none), u2 > 0(none) no shock wave , no rarefaction wave
	#velLeft = 1000.0
	#velRight = -1000.0




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
