import getopt,sys
from constants import problemType
if problemType == "soundwave":
	from model_soundwave import Model
elif problemType == "riemann":
	from model_riemann import Model
else:
	print("problemtype %s not implemented" % problemType)
	sys.exit(0)


def usage():
	print("Usage: python main.py [--timeEnd=<timeEnd>]")



try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["help",  "timeEnd="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
timeEnd = None
for o, a in opts:
	if o in("--timeEnd"):
		timeEnd = float(a)
	elif o in ("--help"):
		usage()
		sys.exit()
	else:
		print("option %s not recognized " % o)
		usage()
if timeEnd is None:
	from constants import timeEnd
m = Model()
m.mainLoop(timeEnd)
import time
time.sleep(5)


