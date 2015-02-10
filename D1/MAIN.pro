
; pro MAIN


; ---------------------------------------------------------
;   -> This is the decision and task-distribution core of the code
;   -> It carries out a *limited* number of assignments itself. 
; ---------------------------------------------------------


; -------------
; COMMON BLOCKS 
; -------------

@params_common.pro
@grid_common
@init_common

; initialize some of the variables in the common blocks
gam=5d0/3d0


; -------------------------------
; READ PARAMETERS FROM INPUT FILE
; -------------------------------

; OPEN THE INPUT DATA FILE
filename='data.dat'
openr,unit,/get_lun,filename
strdum='' 

; READ IN THE PARAMETERS FOR THE NUMERICAL MESH and CREATE THE MESH
readf,unit,nint  & npz=nint+2        
readf,unit,z0,zf  
zz = dindgen(npz)/nint * (zf-z0) + z0  
zz = zz - (zz[1]-zz[0])/2d0

; READ IN MAX ITERATIONS, MAX TIME, OUTPUT PERIODICITY
readf, unit, itmax     
readf, unit, timeend 
readf, unit, plotinterv, storeinterv 

readf, unit, strdum ; dummy reading

; READ IN THE PARAMETERS OF THE INITIAL CONDITION
readf, unit,form='(a25)',strdum & inittype=strtrim(strdum,2)
readf, unit,form='(a25)',strdum &    shape=strtrim(strdum,2)  
readf, unit, amp    

; READ IN THE CFL PARAMETER
readf, unit, fcfl

free_lun,unit ; close and deallocate the unit

; --------------------------------------
; READ PARAMETERS FROM INPUT FILE -- END
; --------------------------------------
firsttime=1


; -------------------------
;  INITIAL CONDITIONS 
; -------------------------

initcond

  momz0 = rho0*vz0 
 energ0 = pres0/(gam-1d0) + rho0*vz0*vz0/2d0 

; do you want to plot the initial conditions for checking? 
; if so: call some plotting routine here. 

; -------------------------
;  INITIAL CONDITIONS - end
; -------------------------


 
; ------------------------------------------------------------
; DECLARE / INITIALIZE VARIABLES JUST BEFORE STARTING THE LOOP
; ------------------------------------------------------------

     vz = vz0   &  rho    = rho0    & pres = pres0  
   momz = momz0 &  energ  = energ0
  rhon  = rho   &  momzn  = momz    & energn = energ
  rhonn = rho   &  momznn = momz    & energnn = energ

it=0L & time=0d0 

; --------------------------------------------------------
; INITIALIZE VARIABLES JUST BEFORE STARTING THE LOOP - end
; --------------------------------------------------------


; ---------------
; BIG LOOP BEGINS
; ---------------

while 1 do begin

;   CALCULATE FLUXES FROM DENSITIES
    fluxes,rho,momz, energ, mflz, momflzz, energflz

;   TIMESTEP
    dt=timestep(rho,momz,pres)

;   UPDATE ACROSS TIMESTEP

    update, dt, rho,momz,energ, mflz,momflzz,energflz, rhonn,momznn,energnn

;   BOUNDARY CONDITIONS
    bcs,'periodic','left',rhonn     & bcs,'periodic','right',rhonn
    bcs,'periodic','left',momznn    & bcs,'periodic','right',momznn
    bcs,'periodic','left',energnn   & bcs,'periodic','right',energnn       

;   CALCULATE PRIMITIVE VARIABLES FROM THE DENSITIES    
    dens2prim, rhonn,momznn,energnn, vznn, presnn

;   EXCHANGE NEW AND OLD VARIABLES
    rho = rhonn & momz = momznn & energ = energnn & pres = presnn  & vz = vznn 

    it = it+1         ; alternative (like in C): it++
    time = time + dt  ; alternative (like in C): time += dt

;   STORE RESULTS
;   write here some command if you want to store results in a file
    if it mod storeinterv eq 0 then begin
        print,'hello'
    endif

;   PLOT RESULTS 
    if it mod plotinterv eq 0 or it eq 0 or it eq itmax then begin
        draw,rho,pres,vz, time, it
    endif

;   CHECK IF THE CALCULATION MUST BE FINISHED
    if it ge itmax or time ge timeend then break
endwhile

; --------------
; BIG LOOP - end
; --------------

print,'program finished.'+$
      '     it= '+ strtrim(string(it),2) +$
      '; time = '+ strtrim(string(time),2)

end
