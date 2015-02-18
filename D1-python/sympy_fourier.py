from sympy import exp , fourier_transform, symbols, exp, pi, cos


z, zc, k0, W , z0, zf, k = symbols("z zc k0 W z0 zf k")

print(fourier_transform(exp(-(z - zc)**2 / W**2) *  cos( 2 * pi * k0 * (z-z0) / (zf - z0)) , z, k))
