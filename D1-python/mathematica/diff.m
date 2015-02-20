
$Assumptions = {Element[{k0,z0,zf,zc,W}, Reals], k0>0, z0>0, zf >0 , zf >0, zc >0, W>0 }

f1[m_]:=((E^(((4*I)*k0*Pi*z0*(zc + zf))/(z0 - zf)^2) + E^((4*k0*Pi*(m*Pi*W^2 + I*(z0^2 + zc*zf)))/(z0 - zf)^2))*Sqrt[Pi]*W)/
  (2*E^((Pi*(k0^2*Pi*W^2 + m*(m*Pi*W^2 - (2*I)*zc*(z0 - zf)) + 2*k0*(m*Pi*W^2 + I*(z0 + zc)*(z0 + zf))))/(z0 - zf)^2))
f2[k_]:=((E^(-((Pi*(k0^2*Pi*W^2 - 2*k0*(k*Pi*W^2 - I*z0 + I*zc)*(z0 - zf) + k*(k*Pi*W^2 + (2*I)*zc)*(z0 - zf)^2))/(z0 - zf)^2)) + 
      E^(-((Pi*(k0^2*Pi*W^2 + 2*k0*(k*Pi*W^2 - I*z0 + I*zc)*(z0 - zf) + k*(k*Pi*W^2 + (2*I)*zc)*(z0 - zf)^2))/(z0 - zf)^2)))*Sqrt[Pi]*W)/2

Print[FullSimplify[f1[k]-f2[k/(zf-z0)]]]
