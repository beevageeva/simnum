#empiric determination of rwPoint and shPoint


def getFirstIndexDifferentLeft(y, delta, default=None):
	#I know they have it
	#if not hasattr(y, '__len__') or len(y) == 0:
	#	return None
	leftVal = y[0]
	for i in range(1, len(y)):
		#because of the smooth functions(for initial conditions: see initcond_riemann.py) I can't test for leftVal == y[i]
		if abs(leftVal - y[i])>delta * leftVal:
			return i
	return default


def getFirstIndexDifferentRight(y, delta, default=None):
	#I know they have it
	#if not hasattr(y, '__len__') or len(y) == 0:
	#	return None
	rightVal = y[len(y) - 1]
	for i in range(len(y) - 1, -1, -1):
		#print("abs dif RIGH %e rightVal = %e , delta * rightVal = %e"  % (abs(rightVal - y[i]), rightVal, delta * rightVal ))	
		if abs(rightVal - y[i])>delta * rightVal:
			return i
	return len(y) - 1 if default is None else default


#def getIndexRightAlmost0(y, delta, startIndex = 0):
#	for i in range(startIndex, len(y)):
#		#print("i = %d, abs(y[i])" % i)
#		#print(abs(y[i]))
#		#print("delta")
#		#print(delta)
#		if abs(y[i])<delta :
#			return i
#	return 0


def getFirstIndexConstant(y, startIndex, delta):
	for i in range(startIndex, len(y)-1):
		#always decreasing???? TODO
		if(y[i-1] - y[i]< delta * 0.05 * y[startIndex]):
			print("startIndex = %d, returning %d" % (startIndex, i+1))
			return i+1
			



