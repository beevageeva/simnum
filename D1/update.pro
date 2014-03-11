
pro update,dt,  rho,momz,energ,     mflz,momflzz,energflz, $
                rhonn,momznn,energnn


; ----------------------------------------------------------------------
; ROUTINE UPDATE 
;
;    PURPOSE:  Calculate variables at timestep n+1
;
;    INPUT ARGUMENTS:  - the densities at timestep n, namely rho, momz, energ
;                      - the fluxes mflz, momflzz, energflz
;
;    COMMON BLOCKS: Note that some necessary input is passed via common
;                       blocks (like the grid parameters and array zz, etc)
;
;    OUTPUT:  the densities at timestep n+1, namely rhonn, momznn, energnn
; ----------------------------------------------------------------------




@params_common.pro
@grid_common


; CARRY OUT HERE THE OPERATIONS NEEDED TO OBTAIN THE DENSITIES rhonn, momznn
; and energnn corresponding to the n+1 timestep

; IMPORTANT NOTE: no boundary conditions should be applied here, just the
; input values for the densities and fluxes

end
