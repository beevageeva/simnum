from sound_wave_params import periodic
if periodic == "x":
	from constants import z0_0 as z0
	from constants import zf_0 as zf
if periodic == "y":
	from constants import z0_1 as z0
	from constants import zf_1 as zf
else:  
	#by default import from dim0
	from constants import z0_0 as z0
	from constants import zf_0 as zf

from math import pi
wl = zf - z0
phi = pi / 6.0 - 2.0 * pi * z0 / (zf - z0)
