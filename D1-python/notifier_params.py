notifierType = "visual"
nstepsPlot = 10 
#nstepsPlot = 1

#The following are used for the notifier(taken from notifier_params as they were wavesound specific)
plotPresCurve = False
#plotPresCurve = True
#plotRhoCurve = False
plotRhoCurve = True
plotVelCurve = False
#plotPresAn = False
#plotRhoAn = False
#plotVelAn = False
plotPresAn = True
plotRhoAn = True
plotVelAn = True
markPoints = False
#markPoints = True
plotVelFFT = True
#plotVelFFT = False
#plotVelFFTAnal=False
plotVelFFTAnal=True
plotPresFFT = True
#plotPresFFT = False

from soundwave_boundary_conditions import periodicType 
if periodicType == "refl": 
	plotPresAn = False
	plotRhoAn = False
	plotVelAn = False

from soundwave_perturbation_params import perturbationType
if perturbationType != "one":
	plotVelFFTAnal = False
else:	
	from soundwave_perturbation_params import functiontype
	from soundwave_medium_params import mediumType
	if functiontype!="wavepacket" or mediumType!="homog":
		plotVelFFTAnal = False
 
