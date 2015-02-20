
(*
Print[Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity}, Assumptions -> {Element[{k0,z0,zf,zc,W}, Reals], k0>0, z0>0, zf >0 , zf >0, zc >0, W>0   }]]
*)
$Assumptions = {Element[{k0,z0,zf,zc,W}, Reals], k0>0, z0>0, zf >0 , zf >0, zc >0, W>0 }


(* export is a crap!
Print[ExportString[ Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ], "TeX" ]]

Print[ExportString[FullSimplify[ Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]], "TeX" ]]
*)

(* with exp 
*)
Print[Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]]
Print[Integrate[Exp[- (z-zc)^2 / W^2] Exp[-2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]]

(*
fullSimplfy makes no difference

Print[FullSimplify[ Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]]]

*)

(* for python! 
Print[ComplexExpand[Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]]]
*)


(* with cos Print[Integrate[Exp[- (z-zc)^2 / W^2] Cos[2 Pi  k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]]

Print[FullSimplify[Integrate[Exp[- (z-zc)^2 / W^2] Cos[2 Pi  k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity} ]]]

*)

