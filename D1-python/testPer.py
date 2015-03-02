from constants import z0, zf

def getPeriodicXi2(xval, a, b):
	p = float(b - a)
	k = int((xval-a)/p)
	print("res1=%e" % xval)
	print("k=%d" %k)
	if k%2==1:
		res = 2.0*b - xval + k*p
	else:
		res = xval - k * p
	print("res2=%e" % res)
	if(res < a):
		res+=p
	if(res > b):
		res-=p
	print("res3=%2.2f, k=%d" % (res,k))	
	return res


def getPeriodicX(xval, a=z0, b=zf):
	if xval > b:
		return 2.0*b - xval
	if xval < a:
		return 2.0*a -xval
	return xval


print("z0=%e, zf=%e" % (z0,zf))
print(getPeriodicX(2, z0 ,zf))
