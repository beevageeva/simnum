
(*
Print[Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity}, Assumptions -> {Element[{k0,z0,zf,zc,W}, Reals], k0>0, z0>0, zf >0 , zf >0, zc >0, W>0   }]]
*)


Print[ExportString[ Integrate[Exp[- (z-zc)^2 / W^2] Exp[2 Pi I k0 (z - z0)/ (zf - z0)] Exp[- 2 Pi I m z / (zf - z0)],  {z, -Infinity, Infinity}, Assumptions -> {Element[{k0,z0,zf,zc,W}, Reals], k0>0, z0>0, zf >0 , zf >0, zc >0, W>0   }], "TeX" ]]

