gamma = 5.0/3

#problemType="soundwave" #problemType may be soundwave or riemann
problemType = "soundwave"
#problemType = "riemann"
if problemType == "soundwave":
	z0 = 3.1
	zf = 7.4
elif problemType == "riemann":
	z0 = -5
	zf = 10
	from riemann_params import riemann_problemType
	if riemann_problemType == "complete":
		gamma = 1.4

#nint = 64
#nint=128
#nint =  256
nint = 1024
#nint = 32

timeEnd = 1.0 #if not set as program argument it's taken from here

#schemeType = "lf"  # scheme type may be lf(Lax - Fr) or fg (first generation)
schemeType = "fg"  
if schemeType == "lf":
	fcfl = 0.99 #use this for lax - fr scheme type
elif schemeType == "fg":
	fcfl = 0.97#use this for first generation scheme
	bcStep = "interm" 
	#bcStep = "final" 


