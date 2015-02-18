from sympy import exp , fourier_transform, symbols, exp, pi, cos, integrate


z, zc, k0, W , z0, zf, k, n = symbols("z zc k0 W z0 zf k n")

print(integrate( cos(2 * pi * n * z / (zf - z0)) * exp(-(z - zc)**2 / W**2) *  cos( 2 * pi * k0 * (z-z0) / (zf - z0)) , z) )



