
(*
Print[Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity}, Assumptions -> {Element[{k0,z0,zf,zc,W}, Reals], k0>0, z0>0, zf >0 , zf >0, zc >0, W>0   }]]
*)
$Assumptions = {Element[{k0,z0,zf,zc,W}, Reals], k0>0, z0>0, zf >0 , zf >0, zc >0, W>0 }


(* export is a crap!
Print[ExportString[ Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ], "TeX" ]]

Print[ExportString[FullSimplify[ Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]], "TeX" ]]
*)


Print[WolframAlpha["Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]", IncludePods -> All, PodStates -> {"Step-by-step solution", "Show all steps"}]]

(*  fullsimplify will output the same!
Print[FullSimplify[ Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]]]
*)
