
PRO INITCOND

; ----------------------------------------------------------------------
; ROUTINE INITCOND
;
;    PURPOSE:  Calculate the arrays of density, velocity and pressure 
;                 at time t=0. 
;    INPUT ARGUMENTS: they are passed via common block init_common
;              inittype:   char variable with the choice of initial condition
;                 shape:   char variable with some subsidiary choice
;                 ** please fill in the list with your choices **
;            
;    COMMON BLOCKS: Note that some additional necessary input is passed via
;                   other common blocks (like the grid parameters,  etc)
;
;    OUTPUT:  the arrays rho0, vz0, pres0 (common block)
; ----------------------------------------------------------------------


@params_common
@grid_common
@init_common

if inittype eq 'sound wave' then begin

    if shape eq 'sine' then begin
     ; *** YOU HAVE TO WRITE HERE THE PROPER DEFINITIONS ***
        rho0 = 1d0 + 0d0*zz
         vz0 =       0d0*zz
       pres0 = 1d0 + 0d0*zz
   endif

    if shape eq 'some_other_shape' then begin
        print,'routine initcond: the some_other_shape shape is not defined yet'
        stop
    endif

endif

if inittype eq 'other' then begin
    print,'routine initcond: no other condition is implemented yet'
    stop
endif

end
