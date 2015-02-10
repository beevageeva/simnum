from prac2 import getPeriodicX, getPeriodicX2, getPeriodicXInt, getPeriodicX2Int
import numpy as np



#c = 1.2
#x = np.linspace(-c,c,65)
#difx = x - 25.88
#res = []
#res2 = []	
#res3 = []	
#res4 = []	
#
#for xv in difx:
#	res.append(getPeriodicX(xv, c))
#	res2.append(getPeriodicXInt(xv, -c, c))
#	res3.append(getPeriodicX2(xv, c))
#	res4.append(getPeriodicX2Int(xv, -c, c))
#
#print("res1 = ")
#print(res)
#print("res2 = ")
#print(res2)
#print("res3 = ")
#print(res3)
#print("res4 = ")
#print(res4)


a = 2
b = 3
x = np.linspace(a,b,65)
difx = x - 28.155
res = []
res2 = []	

for xv in difx:
	res.append(getPeriodicXInt(xv, a, b))
	res2.append(getPeriodicX2Int(xv, a, b))

print("res1 = ")
print(res)
print("res2 = ")
print(res2)
