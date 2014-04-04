from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(1,2,1, projection='3d')
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
print(surf)
print(dir(surf))
ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

ax2 = fig.add_subplot(1,2,2, projection='3d')
surf2 = ax2.plot_surface(X, Y, Z*10, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
plt.draw()
plt.show(block=False)

import time
time.sleep(5)
#TODO
#surf.set_array(np.multiply(2.0, Z).flatten())
#this will go to updateValues in visual_plot
ax.cla()
Z = np.multiply(2.0, Z)
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)


plt.draw()
print("changed")
time.sleep(5)
#plt.show()


