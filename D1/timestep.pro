
FUNCTION TIMESTEP,rho,momz,pres


; ----------------------------------------------------------------------
;
; FUNCTION TIMESTEP
;
;    PURPOSE: Calculate the timestep to guarantee numerical stability
;
;    INPUT ARGUMENTS: those necessary to calculate the sound speed, namely
;           rho, momz, pres
; 
;    COMMON BLOCKS: Note that some necessary input is passed via common
;                       blocks (like gam or the grid zz)
;    OUTPUT:  the delta t 
; 
; ----------------------------------------------------------------------



@params_common.pro
@grid_common

cs = sqrt(gam*pres/rho)
vcharac1 = abs(momz/rho + cs)
vcharac2 = abs(momz/rho - cs)

dt = fcfl* (zz(1)-zz(0)) / max([vcharac1, vcharac2])

return, dt

end
