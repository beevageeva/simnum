gamma = 5.0/3
#gamma = 1.4 #this is for complete riemann problem example
#problemType="soundwave" #problemType may be soundwave or riemann
problemType = "riemann"
if problemType == "soundwave":
	z0 = 3.1
	zf = 7.4
elif problemType == "riemann":
	z0 = -5
	zf = 10
#nint = 64
#nint = 128
nint = 1024
timeEnd = 1.0 #if not set as program argument it's taken from here
verbose = False
#schemeType = "lf"  # scheme type may be lf(Lax - Fr) or fg (first generation)
schemeType = "fg"  
#fcfl = 0.99 #use this for f-l scheme type
fcfl = 0.97#use this for first generation scheme
