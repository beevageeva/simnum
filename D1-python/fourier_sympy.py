from sympy import symbols, fourier_transform, exp, cos

z, km, W, zc, z0, k = symbols('z km W zc z0 k')

print(fourier_transform( exp(-(z-zc)**2 / W**2) * cos(km * (z-z0)), z, k  ))
