perturbationType="one"
#perturbationType="superposition"

if perturbationType == "one":
	A = 3.0 * 10.0 ** (-4)
	#A = 5.0 * 10.0 ** (-2)
	#A = 10.0 ** (-2) #this will work fine as well with fg scheme
	
	#functiontype = 'sine'
	#functiontype = 'gauss'
	functiontype = 'wavepacket'
	#functiontype = 'defined'


elif perturbationType == "superposition":
	A = 3.0 * 10.0 ** (-4)
	#init_functions_generation = [{'csSign':1, 'A': A, 'functiontype': 'sine'}, {'csSign':-1, 'A': 0.5*A, 'functiontype': 'sine'}] #SUPERPOSITION wave travelling right with amp A and left with amp 0.5 * A
	#init_functions_generation = [{'csSign':-1, 'A': A, 'functiontype': 'sine'}, {'csSign':1, 'A': 0.5*A, 'functiontype': 'sine'}] #SUPERPOSITION wave travelling left with amp A and right with amp 0.5 * A
	#init_functions_generation = [{'csSign':-1, 'A': A, 'functiontype': 'sine'}, {'csSign':1, 'A': A, 'functiontype': 'sine'}] #SUPERPOSITION wave travelling left with amp A and right with amp  A
	#init_functions_generation = [{'csSign':-1, 'A': A, 'functiontype': 'sine'}] #travelling left
	#init_functions_generation = [{'csSign':1, 'A': A, 'functiontype': 'sine'}] #travelling right
	#init_functions_generation = [{'csSign':-1, 'A': A, 'functiontype': 'gauss'}, {'csSign':1, 'A': A, 'functiontype': 'gauss'}] #SUPERPOSITION wave travelling left with amp A and right with amp  A but gauss
	#init_functions_generation = [{'csSign':-1, 'A': A, 'functiontype': 'wavepacket'}, {'csSign':1, 'A': A, 'functiontype': 'wavepacket'}] #SUPERPOSITION wave travelling left with amp A and right with amp  A but gauss
	init_functions_generation = [{'csSign':-1, 'A': A, 'functiontype': 'wavepacket'}, {'csSign':1, 'A': A, 'functiontype': 'gauss'}] #SUPERPOSITION wave travelling left with amp A and right with amp  A but gauss

