
PRO BCS,type,site,varr

; ----------------------------------------------------------------------
; ROUTINE BCS
;
;    PURPOSE:  Calculate the boundary conditions of a variable
;    INPUT ARGUMENTS:  PLEASE, FILL THIS IN
;
;    COMMON BLOCKS: Note that some necessary input is passed via common
;                       blocks (like the grid parameters and array zz, etc)
;
;    OUTPUT:  PLEASE, FILL THIS IN
; ----------------------------------------------------------------------


@grid_common


; INITIAL CHECK: type of boundary condition and site
;
bcscond = type ne 'periodic' or (site ne 'left' and site ne 'right')
if bcscond then print,'bcs.pro: bc type or site or not included'


; the if-clause for type is introduced for possible future extensions
if type eq 'periodic' then if site eq 'right' then varr(npz-1) = varr(1)
if type eq 'periodic' then if site eq 'left'  then varr(0) = varr(npz-2)

end

