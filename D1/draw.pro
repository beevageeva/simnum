pro draw,rho,pres,vz,time,it

@params_common.pro
@grid_common
@init_common

if firsttime then begin
    print,'this is a trick to plot something or do some other operation'
    print,'only the first time the program enters this routine'
    firsttime = 0
endif

plot,[0,1],[0,0],yr=[0,1]
xyouts,/norm,0.3,0.6,'PLEASE ... '
xyouts,/norm,0.3,0.5,'BUILD A draw.pro ROUTINE'
xyouts,/norm,0.3,0.4,'THAT PLOTS THE SOLUTIONS YOU OBTAIN'
xyouts,/norm,0.8,0.8,'it='+strtrim(string(it),2)
wait,0.1


end



