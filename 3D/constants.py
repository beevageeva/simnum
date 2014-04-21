gamma = 5.0/3

#problemType="soundwave" #problemType may be soundwave or riemann
problemType = "soundwave"
if problemType == "soundwave":
	z0_0 = 1.2
	zf_0 = 4.3
	z0_1 = 1.2
	zf_1 = 4.3
	z0_2 = 1.2
	zf_2 = 4.3
#	z0_0 = 3.1
#	zf_0 = 7.4
#	z0_1 = 3.1
#	zf_1 = 7.4
	from sound_wave_params import argType
	if(argType == "r"):	
		#in case of periodicity "r" change interval in order to have the whole interval for periodicity in the diagonal direction
		#put domain with center in 0 and make it square (min)
		lhalf = 0.5 *  min((zf_1 - z0_1),(zf_0 - z0_0) )
		z0_0 = z0_1 = -lhalf
		zf_0 = zf_1 = lhalf
	
elif problemType == "riemann":
	z0_0 = z0_1 = -5
	zf_0 = zf_1 = 10
	from riemann_params import riemann_problemType
	if riemann_problemType == "complete":
		gamma = 1.4

nint = 64
#nint =  256
#nint = 1024
#nint = 32

timeEnd = 1.0 #if not set as program argument it's taken from here

verbose = False

#schemeType = "lf"  # scheme type may be lf(Lax - Fr) or fg (first generation)
schemeType = "fg"  
if schemeType == "lf":
	fcfl = 0.99 #use this for lax - fr scheme type
elif schemeType == "fg":
	fcfl = 0.97#use this for first generation scheme