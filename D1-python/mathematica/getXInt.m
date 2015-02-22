(* <<JavaGraphics` *)



(* dens[z_, rho00_, rho01_, ze_, we_]:= rho00 + 0.5 * (rho01-rho00) * (1 + Tanh[(z-ze)/we]) *)
dens[z_]:= rho00 + 0.5 * (rho01-rho00) * (1 + Tanh[(z-ze)/we])
intX[z_]:= FullSimplify[Integrate[Sqrt[dens[z, rho00, rh01, ze, we]], z]]

Print[intX[z]]

(*
Print[FullSimplify[Integrate[Sqrt[dens[z, rho00, rh01, ze, we]], z]]]
intX[z_, rho00_, rho01_, ze_, we_]:= FullSimplify[Integrate[Sqrt[dens[z, rho00, rh01, ze, we]], z]]
*)

(*
intX[z_, rho00_, rho01_, ze_, we_]:=Sqrt[rh01]*we*ArcTanh[(0.7071067811865475*Sqrt[rh01 + rho00 + (rh01 - 1.*rho00)*Tanh[(z - 1.*ze)/we]])/Sqrt[rh01]] -  Sqrt[rho00]*we*ArcTanh[(0.7071067811865475*Sqrt[rh01 + rho00 + (rh01 - 1.*rho00)*Tanh[(z - 1.*ze)/we]])/Sqrt[rho00]]

(* Print[Solve[gp intX[z, rho00, rho01, ze, we] == t + C1, z]] *)
Print[Reduce[gp intX[z, rho00, rho01, ze, we] == t + C1, z]]
i*)
