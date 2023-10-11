

pro carrington_maps

  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

  dir  ='/media/Data1/data1/DATA/ldem_files/'

  rmin = 1.0
  rmax = 1.3
  r0A  = [1.025,1.065,1.105]
  
  

; 3 BANDS
  file ='LDEM._CR2082_euvi.A_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' & suffix3B = 'CR2082_HOLLOW_3Bands'
  read_ldem,file,/ldem,/gauss1,dir=dir
  ratio = sfbe/tfbe
  R     = total( abs(1.-ratio), 4 ) / float(nband)
  ZDA3b   = where(demc eq -999.)
  CNS3b   = where(demc ne -999.  AND R gt 0.25)
  Ne3b = N_e
  Tm3b = Tm
  R3b  = R

  
; 2 BANDS
  file ='LDEM._CR2082_euvi.A_Hollow_2Bands_gauss1_lin_Norm-median_singlStart' & suffix2B = 'CR2082_HOLLOW_2Bands'
  read_ldem,file,/ldem,/gauss1,dir=dir
  ratio = sfbe/tfbe
  R     = total( abs(1.-ratio), 4 ) / float(nband)
  ZDA2b   = where(demc eq -999.)
  CNS2b   = where(demc ne -999.  AND R gt 0.25)
  Ne2b = N_e
  Tm2b = Tm
  R2b  = R

  
  find_min_max,Ne3b,r0A,rmin,rmax,nr,minS_Ne3b,maxS_Ne3b
  find_min_max,Ne2b,r0A,rmin,rmax,nr,minS_Ne2b,maxS_Ne2b
  find_min_max,Tm3b,r0A,rmin,rmax,nr,minS_Tm3b,maxS_Tm3b
  find_min_max,Tm2b,r0A,rmin,rmax,nr,minS_Tm2b,maxS_Tm2b
  minS_Ne = r0A * 0.
  maxS_Ne = r0A * 0.
  minS_Tm = r0A * 0.
  maxS_Tm = r0A * 0.
  for i = 0, n_elements(r0A) - 1 do begin
     minS_Ne(i) = min ( [minS_Ne3b(i),minS_Ne2b(i)])
     maxS_Ne(i) = max ( [maxS_Ne3b(i),maxS_Ne2b(i)])
     minS_Tm(i) = min ( [minS_Tm3b(i),minS_Tm2b(i)])
     maxS_Tm(i) = max ( [maxS_Tm3b(i),maxS_Tm2b(i)])
  endfor
  
 
   
  superhigh=0.25 & p=where(demc ne -999. AND R2b ge superhigh) & if p(0) ne -1 then R2b(p)=superhigh
  superlow=0.01  & p=where(demc ne -999. AND R2b le superlow)  & if p(0) ne -1 then R2b(p)=superlow
  superhigh=0.25 & p=where(demc ne -999. AND R3b ge superhigh) & if p(0) ne -1 then R3b(p)=superhigh
  superlow=0.01  & p=where(demc ne -999. AND R3b le superlow)  & if p(0) ne -1 then R3b(p)=superlow
  
  
  if ZDA2b(0) ne -1 then begin
     Ne2b(ZDA2b)=0.
     Tm2b(ZDA2b)=0.
     R2b (ZDA2b)=0.
  endif
 if ZDA3b(0) ne -1 then begin
    Ne3b(ZDA3b)=0.
    Tm3b(ZDA3b)=0.
    R3b (ZDA3b)=0.
  endif


  
  if CNS2b(0) ne -1 then begin
     Ne2b(CNS2b)= 1.e16
     Tm2b(CNS2b)= 1.e12
     R2b (CNS2b)= 0.25
  endif
 if CNS3b(0) ne -1 then begin
     Ne3b(CNS3b)= 1.e16
     Tm3b(CNS3b)= 1.e12
     R3b (CNS3b)= 0.25
  endif

 
 xdisplay,map=Ne2b,file='Ne_'+suffix2b,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 34, titulo='Ne [10!U8!Ncm!U-3!N]' ,units=1.e8, minA=minS_Ne/1.e8, maxA=maxS_Ne/1.e8,/add_bw
 xdisplay,map=Ne3b,file='Ne_'+suffix3b,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 34, titulo='Ne [10!U8!Ncm!U-3!N]' ,units=1.e8, minA=minS_Ne/1.e8, maxA=maxS_Ne/1.e8,/add_bw
 xdisplay,map=tm2b,file='Tm_'+suffix2b,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 34, titulo='Tm [MK]'              ,units=1.e6, minA=minS_Tm/1.e6, maxA=maxS_Tm/1.e6,/add_bw
 xdisplay,map=tm3b,file='Tm_'+suffix3b,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 34, titulo='Tm [MK]'              ,units=1.e6, minA=minS_Tm/1.e6, maxA=maxS_Tm/1.e6,/add_bw

 xdisplay,map=R2b, file='R_' +suffix2b,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 12, titulo='R'                               , minA=(r0A*0.+1.E-6), maxA=(r0A*0.+0.25)
 xdisplay,map=R3b, file='R_' +suffix3b,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 12, titulo='R'                               , minA=(r0A*0.+1.E-6), maxA=(r0A*0.+0.25)
  


  return
end
