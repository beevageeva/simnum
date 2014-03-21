#empiric determination of rwPoint and shPoint


def getFirstIndexDifferentLeft(y, delta, default):
	if not hasattr(y, '__len__') or len(y) == 0:
		return None
	leftVal = y[0]
	for i in range(1, len(y)):
		#because of the smooth functions(for initial conditions: see initcond_riemann.py) I can't test for leftVal == y[i]
		if abs(leftVal - y[i])>delta * leftVal:
			return i
	return default


def getFirstIndexDifferentRight(y, delta, default):
	if not hasattr(y, '__len__') or len(y) == 0:
		return None
	rightVal = y[len(y) - 1]
	for i in range(len(y) - 1, -1, -1):
		if abs(rightVal - y[i])>delta * rightVal:
			return i
	return default
