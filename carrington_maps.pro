; Este script sirve para hacer mapas lat-lon de N_e y T_m DEMT.
; Explora una nueva escala de color y satura ZDAs en negro y CNSs en
; blanco.
; Ademas, hace histogramas en una box.

pro carrington_maps,NFS=NFS,write=write

  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

  root_dir = '/data1/'
  if keyword_set(NFS) then root_dir='/media/Data1/data1/'
  dir  =root_dir+'DATA/ldem_files/'
  
  file ='LDEM._CR2082_euvi.A_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' & suffix = 'CR2082_HOLLOW_3Bands' 
; file ='LDEM._CR2082_euvi.A_Hollow_2Bands_gauss1_lin_Norm-median_singlStart' & suffix = 'CR2082_HOLLOW_2Bands' 
; file ='LDEM._CR2082_euvi.A_Hollow_2Bands_Tmax2.5MK_gauss1_lin_Norm-median_singlStart' & suffix = 'CR2082_HOLLOW_2Bands_new'  
; file ='LDEM._CR2082_euvi.A_Hollow_2Bands_0.75-2.0MK_gauss1_lin_Norm-median_singlStart'& suffix = 'CR2082_HOLLOW_2Bands_new'  
  if not keyword_set(write) then read_ldem,file,/ldem,/gauss1,dir=dir  
  if     keyword_set(write) then begin
     read_ldem,file,/ldem,/gauss1,dir=dir,/write,/mark_AEVs
     RETURN
  endif

  ratio = sfbe/tfbe &  R = total( abs(1.-ratio), 4 ) / float(nband)
                 ZDA   = where(demc eq -999.)
  R_th  = 0.25 & CNS   = where(demc ne -999.  AND R gt R_th)

  Nesat = N_e
  Tmsat = Tm
  Rsat  = R

  superhigh=R_th & p=where(demc ne -999. AND R gt superhigh) & if p(0) ne -1 then Rsat(p)=superhigh
  superlow=0.01  & p=where(demc ne -999. AND R le superlow)  & if p(0) ne -1 then Rsat(p)=superlow
  
  
  if ZDA(0) ne -1 then begin
     Nesat(ZDA)=0.
     Tmsat(ZDA)=0.
     Rsat (ZDA)=0.
  endif

  if CNS(0) ne -1 then begin
     Nesat(CNS)= 1.e16
     Tmsat(CNS)= 1.e12
     Rsat (CNS)= R_th
  endif
  
; r0A=    [1.025,1.065,1.105,1.155,1.205,1.245]
; maxA_Ne=[2.0  ,1.75 ,1.5  ,1.25 ,1.0  ,0.5 ]
  r0A=    [1.025,1.105,1.245]
  maxA_Ne=[2.0  ,1.5  ,0.5 ]

  
  xdisplay,map=Nesat,file='Ne_'+suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 25,$
           titulo='Ne [10!U8!Ncm!U-3!N]', units=1.e8, minA=(r0A*0.+1.E-6), maxA=maxA_Ne, /add_bw
  xdisplay,map=tmsat,file='Tm_'+suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 25,$
           titulo='Tm [MK]',units=1.e6, minA=(r0A*0.+1.E-6), maxA=3.0+r0A*0., /add_bw
  xdisplay,map=Rsat, file='R_' +suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 12,$
           titulo='R',minA=(r0A*0.+1.E-6), maxA=(r0A*0.+0.25), /add_bw

  

  STOP

; =============================================================================
; Make the histograms for the heights set in r0A 
  if CNS(0) ne -1 then begin
    Nesat(CNS)= -666.
    Tmsat(CNS)= -666.
  endif


; Histogram in the North CH
;-----------------------------------------------------------------------------
  minA_Ne=[0.5,0.2,0.05]
  maxA_Ne=[1.5,0.6,0.20]
  minA_Te=[0.5,0.5,0.5]
  maxA_Te=[1.5,1.5,1.5]
  r0A=[1.025,1.105,1.245] 
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
; -----------------------------------------------------------------------------

; Histogram in the Streamer Region
  
  minA_Ne=[0.5,0.2,0.05]*2.
  maxA_Ne=[1.5,0.6,0.20]*2.
  minA_Te=[0.5,0.5,0.5]
  maxA_Te=[2.0,2.0,2.0]
  r0A=[1.025,1.105,1.245] 
  for i=0,n_elements(r0A)-1 DO BEGIN
     height        = r0A[i]
     height_string = strmid(string(height),6,5)
     height_string = (STRJOIN(STRSPLIT(height_string, /EXTRACT,'.'), '_'))
     file   = suffix+'_'+height_string+'_Rsun'
     rmin   = 1.0
     rmax   = 1.3
     dr = (rmax-rmin)/nr
     rad_range = [height-dr/2,height+dr/2]
     lat_range = [-50.,+50.] 
     lon_range = [  0.,360.]
     ;yrange=[0.,0.1]
     xhisto,map=Nesat/1.e8,rmin=rmin,rmax=rmax,nr=nr,nt=nth,rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,titulo='N!de!n [10!U8!Ncm!U-3!N]',file='Ne_'+file,sufijo='Str',$
            mini=minA_Ne[i],maxi=maxA_Ne[i];,yrange=yrange
     ;yrange=[0.,0.2]
     xhisto,map=Tmsat/1.e6,rmin=rmin,rmax=rmax,nr=nr,nt=nth,rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,titulo='T!dm!n [MK]',file='Tm_'+file,sufijo='Str',$
            mini=minA_Te[i],maxi=maxA_Te[i];,yrange=yrange
  ENDFOR


  
; =============================================================================




  
  return
end
