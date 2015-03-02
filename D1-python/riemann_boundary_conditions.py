def  lrBoundaryConditions(array, skip=0):
	n = array.shape[0] - 1
	array = np.insert(array, 0, array[0])
	array = np.insert(array, n+2, array[-1])
	return array
