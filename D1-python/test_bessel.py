from scipy.special import jv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0,16,1000) 
y = jv(0,x)
plt.plot(x,y,"r-")
plt.plot()
plt.show()
