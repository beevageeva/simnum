(* <<JavaGraphics` *)

$Assumptions = {Element[{k0,z0,zf,zc,W}, Reals], k0>0, z0>0, zf >0 , zf >0, zc >0, W>0 }


h[z_, k0_, z0_, zf_, zc_, W_]:= Exp[-(z-zc)^2/W^2] Cos[2 Pi k0 (z - z0) / (zf - z0) ]
(*
hh[z_] := h[z, k0, z0, zf, zc, W]
Print[FullSimplify[FourierTransform[hh[z], z, k]]]
Print[FullSimplify[FourierTransform[hh[z], z, k, FourierParameters->{0,-2 Pi}]] ]
*)


Print[FullSimplify[FourierTransform[h[z, k0, z0, zf, zc, W], z, k, FourierParameters->{0,-2 Pi}]] ]

(*
Print[ComplexExpand[FullSimplify[FourierTransform[h[z, k0, z0, zf, zc, W], z, k, FourierParameters->{0,-2 Pi}]] ]]
*)

(*

ft[k_]:= FullSimplify[FourierTransform[hh[z], z, k]]

ft2[k_]:= FullSimplify[FourierTransform[hh[z], z, k, FourierParameters->{0,-2 Pi}]] 

*)

(*
ft[k_, k0_, z0_, zf_, zc_, W_]:=(E^(-(zc^2/W^2) + ((2*I)*k0*Pi*z0)/(z0 - zf))*(E^((W^2*(I*k + (2*zc)/W^2 + ((2*I)*k0*Pi)/(-z0 + zf))^2)/4) + E^((W^2*(I*k + (2*zc)/W^2 + ((2*I)*k0*Pi)/(z0 - zf))^2)/4 + ((4*I)*k0*Pi*z0)/(-z0 + zf))))/(2*Sqrt[2]*Sqrt[W^(-2)])

ft2[k_, k0_, z0_, zf_, zc_, W_]:= (E^(-(zc^2/W^2) + ((2*I)*k0*Pi*z0)/(z0 - zf))*(E^((W^2*((2*zc)/W^2 - (2*I)*Pi*(k + k0/(z0 - zf)))^2)/4) + E^(((4*I)*k0*Pi*z0)/(-z0 + zf) + (W^2*((2*zc)/W^2 - (2*I)*Pi*(k + k0/(-z0 + zf)))^2)/4))*Sqrt[Pi])/(2*Sqrt[W^(-2)])
((E^(-((Pi*(k0^2*Pi*W^2 - 2*k0*(k*Pi*W^2 - I*z0 + I*zc)*(z0 - zf) + k*(k*Pi*W^2 + (2*I)*zc)*(z0 - zf)^2))/(z0 - zf)^2)) + 
      E^(-((Pi*(k0^2*Pi*W^2 + 2*k0*(k*Pi*W^2 - I*z0 + I*zc)*(z0 - zf) + k*(k*Pi*W^2 + (2*I)*zc)*(z0 - zf)^2))/(z0 - zf)^2)))*Sqrt[Pi]*W)/2

z0 = 3.1
zf = 7.4
k0 = 60.0
(* zc = z0 + (3.0/20.0)*(zf - z0) *)(* second exp of inhom and first new *)
zc = z0
W = 0.05 
(* W = 0.5 *)
numPoints = 1024
(*
k = FindDivisions[{0, 2 * Pi* k0},numPoints ]
y1 = Abs[ft[k, k0, z0, zf, zc, W]]
plot1 = Plot[y1, {k, 0, 2*Pi*k0}]
index1 = Position[y1, Max[y1]]
Print["index1"]
Print[index1]
Print["k1"]
Print[k[[index1]]]
Print[k[[176]]]
Export["plot1.png", plot1]
*)
(* k = FindDivisions[{-4/3 k0,  4/3 k0},numPoints ] *)
k = FindDivisions[{-k0,  k0},numPoints ]
(*
Print[ft2[k, k0, z0, zf, zc, W]]
Print["PEAK"]
Print[ft2[k0, k0, z0, zf, zc, W]]
*)
y2 = Abs[ft2[k, k0, z0, zf, zc, W]]
plot2 = Plot[Abs[ft2[k, k0, z0, zf, zc, W]], {k, -k0, k0}]
index1 = Position[y2, Max[y2]]
Print["index2"]
Print[index1]
Print["k1"]
Print[k[[index1]]]

Print[k[[141]]]
Export["plot2n.png", plot2]

*)

(*
Print["*****************FT"]
(* Print[??ft] 
DumpSave["ft.mx", ft] *)

Print[ft[900]]
Export["fd1.txt", ft]
Export["fd2.txt", ft2]

Print["*****************FT2"]
Print[ft2[60]]
Print[ft2[1]]


z = FindDivisions[{z0, zf},numPoints ]
Print[z]
Print["Values"]
Print[hh[z]]
numFourier = Fourier[hh[z]]


plot1 = ListLinePlot[Abs[numFourier], PlotRange->Full]
plot2 = Plot[Abs[ft[k]], {k, 1,1000}, PlotRange->All]

plot3 = ListLinePlot[Abs[ft[k]], PlotRange->Full]

f1[k_]:=FourierTransform[Sin[10^3 t] Exp[-t/10] UnitStep[t], t, k]
plot4 = LogLogPlot[Abs[f1[k]], {k, 1*^-1, 1*^5}, PlotRange -> All]


Export["plot1.png", plot1]
Export["plot2.png", plot2]
Export["plot3.png", plot3]
Export["plot4.png", plot4]

Print[h[z0, k0, z0, zf, zc, W]]
Print[h[zc, k0, z0, zf, zc, W]]
Print[h[zf, k0, z0, zf, zc, W]]
*)
