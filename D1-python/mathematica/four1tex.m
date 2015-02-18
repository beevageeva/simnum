
h[z_, k0_, z0_, zf_, zc_, W_]:= Exp[-(z-zc)^2/W^2] Cos[2 Pi k0 (z - z0) / (zf - z0) ]
Print[ExportString[FullSimplify[FourierTransform[h[z, k0, z0, zf, zc, W], z, k, FourierParameters->{0,-2 Pi}],Element[{k0,z0,zf,zc,W}, Reals] ] , "TeX"]]
