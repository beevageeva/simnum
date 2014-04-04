import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np


X = np.arange(0, 10, 0.1)
Y = np.arange(0, 10, 0.1)
Z = np.arange(0, 10, 0.1)


fig = plt.figure(figsize=(12,6))


ax = fig.add_subplot(1,2,1, projection='3d')
#ax.plot_surface(X, Y, Z, rstride=4, cstride=4, alpha=0.25)
#p = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=False)
p = ax.plot_wireframe(X, Y, Z, rstride=4, cstride=4)
ax.view_init(30, 45)
#ax.view_init(70, 30)


fig.tight_layout()

plt.draw()
plt.show()
