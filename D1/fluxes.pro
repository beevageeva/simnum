
pro fluxes,rho,momz,energ,     mflz,momflzz,energflz

; ----------------------------------------------------------------------
;
; ROUTINE FLUXES
;
;    PURPOSE: Calculate the fluxes f_m, f_c, f_e 
;
;    INPUT ARGUMENTS: the primitive variables density (rho), momentum (momz)
;                     and total energy (energ)
;
;    COMMON BLOCKS: Note that some necessary input is passed via common
;                       blocks (like the grid parameters and array zz, etc)
;
;    OUTPUT:  the fluxes mflz, momflzz, energflz 
; ----------------------------------------------------------------------



  dens2prim, rho,momz,energ,  vz,pres ; this calculates vz,pres

  mflz = momz   ; mass flux
  momflzz = momz*momz/rho + pres   ; momentum flux
  energflz = (energ + pres) * vz   ; energy flux

end
