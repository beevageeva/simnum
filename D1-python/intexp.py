from sympy import Symbol, symbols, tanh, integrate, exp, oo

z = Symbol("z")

print(integrate(exp(-z**2), (z, -oo, oo) ))
