from scipy.special import hankel1
import matplotlib.pyplot as plt
import numpy as np

def deriv(f):
  n = len(f)
  res = np.zeros(n)
  for i in range(1, n-1):
    res[i] = f[i+1] - f[i-1] / (2*dx)
  res[0] = res[1]
  res[n-1] = res[n-2]
  return res

	
nint = 128
a = 5

x = np.linspace(-a,a,nint+2)
index0=nint/2 
dx = x[1] - x[0]

def h10mod(x):
	res = np.zeros(x.shape)
	xm = 1
	index1=int(nint * (xm+a) / (2.0*a) )	
	index_1=int(nint * (a-xm) / (2.0*a) )	
	for i in  range(index1, len(x)):
		res[i] = hankel1(0,x[i])
	for i in  range(0, index_1+1):
		res[i] = -hankel1(0,x[i])
	for i in  range(index_1+1, index1):
		res[i] = res[index1]
	return res


# Hankel function H1_n(x) on the real line for n=0,1,2,3
h0 = lambda x: hankel1(0,x)
h1 = lambda x: hankel1(1,x)
h2 = lambda x: hankel1(2,x)
h3 = lambda x: hankel1(3,x)
vals = h0(x)
vals2 = h10mod(x)	
#plt.plot(x, vals2, color="k")
#plt.plot(x, vals, color="g")
plt.plot(x, -h1(x).real, color="r")
plt.plot(x, deriv(vals).real, color="b")
plt.plot(x, deriv(vals2).real, color="y")
plt.draw()
plt.show()

