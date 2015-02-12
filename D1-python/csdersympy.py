from sympy import Symbol, symbols, tanh, diff, sqrt

rho00, rho01, z, ze, we  = symbols("rho00 rho01 z ze we")

print(diff(1/sqrt(rho00 + 0.5 * (rho01-rho00) * (1 + tanh((z-ze)/we))), z))
