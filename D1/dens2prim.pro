
pro dens2prim,rho,momz,energ,   vz,pres

; ----------------------------------------------------------------------
; ROUTINE dens2prim
;
;    PURPOSE:  Calculate the primitive variables vz and pres from the
;              mass, momentum and energy densities.
;
;    INPUT ARGUMENTS:  PLEASE, FILL THIS IN
;
;    COMMON BLOCKS: Note that some necessary input is passed via common
;                       blocks 
;
;    OUTPUT:  PLEASE, FILL THIS IN
; ----------------------------------------------------------------------



@params_common.pro

vz = momz/rho
pres = (gam-1d0)*(energ-momz*momz/2d0/rho)  

end
