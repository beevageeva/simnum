
dens[z_]:= rho00 + 0.5 * (rho01-rho00) * (1 + Tanh[(z-ze)/we])
intX[z_]:= FullSimplify[Integrate[Sqrt[dens[z]], z]]


Print[intX[z]]

resX[z_] = -0.9999999999999998*Sqrt[rho00]*we*ArcTanh[(0.7071067811865475*Sqrt[rho00 + rho01 + (-1.*rho00 + rho01)*Tanh[(z - 1.*ze)/we]])/Sqrt[rho00]] + 0.9999999999999998*Sqrt[rho01]*we*ArcTanh[(0.7071067811865475*Sqrt[rho00 + rho01 + (-1.*rho00 + rho01)*Tanh[(z - 1.*ze)/we]])/Sqrt[rho01]] 

Print["----------"]
Print[resX[z]]
Print["----------"]

we = 0.4
z0 = 4.1
zf = 7.2
rho00 = 1
rho01 = 0.2
ze = 0.5 * (z0 + zf)


Print[resX[z]]
Print["----------"]

Print[resX[4]]

