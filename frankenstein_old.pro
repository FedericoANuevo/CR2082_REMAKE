

pro FRANKENSTEIN

  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e

; Path to use de x-tools (SolarTom_idl)  
  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

; Suffix to create the graph file-names 
  suffix='frankenstein'
  
; Directory where the DEMT files are stored
  dir  ='/media/Data1/data1/DATA/ldem_files/' ; IAFE
  dir  ='/data1/DATA/ldem_files/'             ; HOME-OFFICE

; DEMT files
  file_3bands ='LDEM._CR2082_euvi.A_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' 
  file_2bands ='LDEM._CR2082_euvi.A_Hollow_2Bands_gauss1_lin_Norm-median_singlStart'  

; Read DEMT / 3 BANDS (EUVI-3)
  read_ldem,file_3bands,/ldem,/gauss1,dir=dir
  R         = total( abs(1.-sfbe/tfbe), 4 ) / float(nband)
  Ne_3bands = N_e
  Tm_3bands = Tm
  R_3bands   = R
  
; FBEs 
  fbe_171 = tfbe(*,*,*,0)
  fbe_195 = tfbe(*,*,*,1)
  fbe_284 = tfbe(*,*,*,2)

  
; Read DEMT / 2 BANDS (EUVI-2)
  read_ldem,file_2bands,/ldem,/gauss1,dir=dir
  R         = total( abs(1.-sfbe/tfbe), 4 ) / float(nband)
  Ne_2bands = N_e
  Tm_2bands = Tm
  R_2bands   = R
  
; DONDE HAY ZDAs en 284 pero NO en las otras dos bandas
  index = where (fbe_284 le 0. and fbe_171 gt 0. and fbe_195 gt 0.)

; REEMPLAZO EL RESULTADO POR EL DE 2 BANDAS
  N_e = Ne_3bands & N_e (index) = Ne_2bands(index)
  Tm  = Tm_3bands & Tm  (index) = Tm_2bands(index)
  R   = R_3bands  & R   (index) = R_2bands (index)


; DONDE hay CNSs para EUVI-3 pero NO hay CNSs para EUVI-2
  index2 = where (demc ne -999. and R_2bands lt 0.25 and R_3bands ge 0.25)

; REEMPLAZO EL RESULTADO POR EL DE 2 BANDAS
  N_e (index2) = Ne_2bands(index2)
  Tm  (index2) = Tm_2bands(index2)
  R   (index2) = R_2bands (index2)

; =============================== GRAPH MODULE =============================================================


; Heighs to plot
  r0A=[1.025,1.065,1.105,1.205]
  rmin=1.0
  rmax=1.3
  
  ZDA   = where(demc eq -999.)
  CNS   = where(demc ne -999.  AND R gt 0.25)


; arrays para salturar en ZDAs y CNSs  
  Nesat = N_e
  Tmsat = Tm
  Rsat  = R

  find_min_max,Nesat,r0A,rmin,rmax,nr,NeminS,NemaxS,NemedS,NemeanS,NestdvS
  find_min_max,Tmsat,r0A,rmin,rmax,nr,TeminS,TemaxS,TemedS,TemeanS,TestdvS
 
; ZDAs en NEGRO 
  if ZDA(0) ne -1 then begin
     Nesat(ZDA)=0.
     Tmsat(ZDA)=0.
     Rsat (ZDA)=0.
  endif

; CNSs en BLANCO
  if CNS(0) ne -1 then begin
 ; AHORA SATURO EN NEGRO COMO ZDAs    
     Nesat(CNS)= 0.;1.e16
     Tmsat(CNS)= 0.;1.e12
     Rsat (CNS)=0.25
  endif

; Para que el mapa de R, donde R << 1, sea verde oscuro y no negro  
  superlow=0.01  & p=where(demc ne -999. AND R le superlow)  & if p(0) ne -1 then Rsat(p)=superlow


; Make the lat-lon maps using xdisplay   
  xdisplay,map=Nesat,file='Ne_'+suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 25, titulo='Ne [10!U8!Ncm!U-3!N]' ,units=1.e8,/add_bw,$
           minA=(NemeanS-NestdvS*1.2)/1.e8,maxA=(NemeanS+NestdvS*1.5)/1.e8;,minA=(r0A*0.+1.E-6), maxA=[2.0,1.5,1.0]
  xdisplay,map=Tmsat,file='Tm_'+suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 25, titulo='Tm [MK]'              ,units=1.e6,/add_bw,$
           minA=(TemeanS-TestdvS*1.2)/1.e6,maxA=(TemeanS+TestdvS*1.5)/1.e6 ;,minA=(r0A*0.+1.E-6), maxA=[2.0,2.0,2.0]
  xdisplay,map=Rsat ,file='R_' +suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 12, titulo='R', minA=(r0A*0.+1.E-6), maxA=(r0A*0.+0.25)



  return
end
