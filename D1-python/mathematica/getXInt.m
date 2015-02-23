(* <<JavaGraphics` *)



(* dens[z_, rho00_, rho01_, ze_, we_]:= rho00 + 0.5 * (rho01-rho00) * (1 + Tanh[(z-ze)/we]) *)
dens[z_]:= rho00 + 0.5 * (rho01-rho00) * (1 + Tanh[(z-ze)/we])
intX[z_]:= FullSimplify[Integrate[Sqrt[dens[z]], z]]


Print[intX[z]]
(* Print[ExportString[intX[z] , "TeX"]] *)


Print["----------"]

we = 0.4
z0 = 4.1
zf = 7.2
rho00 = 1
rho01 = 0.2
ze = 0.5 * (z0 + zf)

resX[z_] = -0.9999999999999998*Sqrt[rho00]*we*ArcTanh[(0.7071067811865475*Sqrt[rho00 + rho01 + (-1.*rho00 + rho01)*Tanh[(z - 1.*ze)/we]])/Sqrt[rho00]] + 0.9999999999999998*Sqrt[rho01]*we*ArcTanh[(0.7071067811865475*Sqrt[rho00 + rho01 + (-1.*rho00 + rho01)*Tanh[(z - 1.*ze)/we]])/Sqrt[rho01]] 

Print[resX[4]]

(*
Print[FullSimplify[Integrate[Sqrt[dens[z, rho00, rh01, ze, we]], z]]]
intX[z_, rho00_, rho01_, ze_, we_]:= FullSimplify[Integrate[Sqrt[dens[z, rho00, rh01, ze, we]], z]]
*)

(*
intX[z_, rho00_, rho01_, ze_, we_]:=Sqrt[rh01]*we*ArcTanh[(0.7071067811865475*Sqrt[rh01 + rho00 + (rh01 - 1.*rho00)*Tanh[(z - 1.*ze)/we]])/Sqrt[rh01]] -  Sqrt[rho00]*we*ArcTanh[(0.7071067811865475*Sqrt[rh01 + rho00 + (rh01 - 1.*rho00)*Tanh[(z - 1.*ze)/we]])/Sqrt[rho00]]

(* Print[Solve[gp intX[z, rho00, rho01, ze, we] == t + C1, z]] *)
Print[Reduce[gp intX[z, rho00, rho01, ze, we] == t + C1, z]]
i*)
