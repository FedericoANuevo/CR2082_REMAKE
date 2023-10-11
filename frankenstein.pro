pro FRANKENSTEIN,NFS=NFS,WRITE=WRITE

  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

  root_dir = '/data1/'
  if keyword_set(NFS) then root_dir='/media/Data1/data1/'
  dir  =root_dir+'DATA/ldem_files/'

  file_3bands ='LDEM._CR2082_euvi.A_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' 
  file_2bands ='LDEM._CR2082_euvi.A_Hollow_2Bands_gauss1_lin_Norm-median_singlStart'  

  read_ldem,file_3bands,/ldem,/gauss1,dir=dir
    R_3bands = total( abs(1.-sfbe/tfbe), 4 ) / float(nband)
   R0_3bands = reform( (abs(1.-sfbe/tfbe))(*,*,*,0) )
   R1_3bands = reform( (abs(1.-sfbe/tfbe))(*,*,*,1) )
   R2_3bands = reform( (abs(1.-sfbe/tfbe))(*,*,*,2) )
  fbe_3bands = tfbe
   Ne_3bands = N_e
   Tm_3bands = Tm
 demc_3bands = demc   
  ZDA_3bands = where(demc eq -999.)

  read_ldem,file_2bands,/ldem,/gauss1,dir=dir
     R_2bands = total( abs(1.-sfbe/tfbe), 4 ) / float(nband)
   fbe_2bands = tfbe
    Ne_2bands = N_e
    Tm_2bands = Tm
 demc_2bands = demc 
  ZDA_2bands = where(demc eq -999.)
   
; Set threshholds to handle CNSs.
  R_th_3b   = 0.25
  R_th_2b   = 0.25
  R_th_comp = 0.25
  
; Initialize compound results to 3Bands
  N_e =   Ne_3bands
  Tm  =   Tm_3bands
  R   =    R_3bands
 demc = demc_3bands

; Assign 2bands results where 3Bands has ZDA only due to third band (so that 2bands has no ZDA for sure).
  index = where (fbe_3bands(*,*,*,0) gt 0. AND fbe_3bands(*,*,*,1) gt 0. AND fbe_3bands(*,*,*,2) le 0.)
  N_e (index) =    Ne_2bands(index)
  Tm  (index) =    Tm_2bands(index)
  R   (index) =     R_2bands(index)
  demc(index) =  demc_2bands(index)


; Assign 2bands results where 3Bands has CNS and 2Bands do not.
  suffix = 'CR2082_HOLLOW_compound1'
  index = where (demc_3bands ne -999. AND R_3bands gt R_th_3b AND R_2bands lt R_th_2b)
  N_e (index) =    Ne_2bands(index)
  Tm  (index) =    Tm_2bands(index)
  R   (index) =     R_2bands(index)
  demc(index) =  demc_2bands(index)


; Determine ZDA and CNS indexes for compound result
  ZDA   = where(demc eq -999.)
  CNS   = where(demc ne -999.  AND R gt R_th_comp)


  IF keyword_set(write) then begin
     N_e(CNS) = -666.
     Tm (CNS) = -666.
   ; Save the results  
     openw,1,dir+'Ne_'+suffix+'.dat'
     writeu,1,float(N_e)
     close,1
     openw,1,dir+'Te_'+suffix+'.dat'
     writeu,1,float(Tm) 
     close,1
     RETURN
  ENDIF
  
  
; Display compound results
  Nesat = N_e
  Tmsat = Tm
  Rsat  = R

  superhigh = R_th_comp & p=where(demc ne -999. AND R gt superhigh) & if p(0) ne -1 then Rsat(p)=superhigh
  superlow  = 0.01      & p=where(demc ne -999. AND R le superlow)  & if p(0) ne -1 then Rsat(p)=superlow
  
  
  if ZDA(0) ne -1 then begin
     Nesat(ZDA)=0.
     Tmsat(ZDA)=0.
     Rsat (ZDA)=0.
  endif

  if CNS(0) ne -1 then begin
     Nesat(CNS)= 1.e16
     Tmsat(CNS)= 1.e12
     Rsat (CNS)= R_th_comp
  endif


; r0A=    [1.025,1.065,1.105,1.155,1.205,1.245]
; maxA_Ne=[2.0  ,1.75 ,1.5  ,1.25 ,1.0  ,0.5 ]
  r0A=    [1.025,1.105,1.245]
  maxA_Ne=[2.0  ,1.5  ,0.5 ]

  xdisplay,map=Nesat,file='Ne_'+suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 25, titulo='Ne [10!U8!Ncm!U-3!N]' ,units=1.e8, minA=(r0A*0.+1.E-6), maxA=maxA_Ne, /add_bw
  xdisplay,map=tmsat,file='Tm_'+suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 25, titulo='Tm [MK]'              ,units=1.e6, minA=(r0A*0.+1.E-6), maxA=3.0+r0A*0., /add_bw
  xdisplay,map=Rsat, file='R_' +suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 12, titulo='R'                               , minA=(r0A*0.+1.E-6), maxA=(r0A*0.+0.25)

  STOP


; =============================================================================
; Make the histograms for the heights set in r0A 
  if CNS(0) ne -1 then begin
    Nesat(CNS)= -6.
    Tmsat(CNS)= -6.
  endif

  r0A=[1.025,1.105,1.245]
  minA_Ne=[0.5,0.2,0.05]
  maxA_Ne=[1.5,0.6,0.20]
  minA_Te=[0.5,0.5,0.5]
  maxA_Te=[1.5,1.5,1.5]
 
  
  for i=0,n_elements(r0A)-1 DO BEGIN
     height        = r0A[i]
     height_string = strmid(string(height),6,5)
     height_string = (STRJOIN(STRSPLIT(height_string, /EXTRACT,'.'), '_'))
     file   = suffix+'_'+height_string+'_Rsun'
     rmin   = 1.0
     rmax   = 1.3
     dr = (rmax-rmin)/nr
     rad_range = [height-dr/2,height+dr/2]
     lat_range = [+60.,+90.] 
     lon_range = [  0.,360.]
     yrange=[0.,0.1]
     xhisto,map=Nesat/1.e8,rmin=rmin,rmax=rmax,nr=nr,nt=nth,rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,titulo='N!de!n [10!U8!Ncm!U-3!N]',file='Ne_'+file,sufijo='NCH',$
            mini=minA_Ne[i],maxi=maxA_Ne[i],yrange=yrange
     yrange=[0.,0.2]
     xhisto,map=Tmsat/1.e6,rmin=rmin,rmax=rmax,nr=nr,nt=nth,rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,titulo='T!dm!n [MK]',file='Tm_'+file,sufijo='NCH',$
            mini=minA_Te[i],maxi=maxA_Te[i],yrange=yrange
  ENDFOR
; =============================================================================

  
  
  return
end
