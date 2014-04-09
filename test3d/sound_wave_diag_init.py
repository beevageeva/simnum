from math import pi

def getPeriodicKx():
	from constants import z0_0, zf_0
	return 	2 * pi / (zf_0 - z0_0)

def getPeriodicKy():
	from constants import z0_1, zf_1
	return 	2 * pi / (zf_1 - z0_1)

