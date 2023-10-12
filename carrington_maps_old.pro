
; Este codigo usa la escala de color historica (de los papers) y
; satura en negro tanto ZDAs como CNSs

pro carrington_maps

  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

  dir  ='/media/Data1/data1/DATA/ldem_files/'
  file ='LDEM._CR2082_euvi.A_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' & suffix = 'CR2082_HOLLOW_3Bands' 
  ;file ='LDEM._CR2082_euvi.A_Hollow_2Bands_gauss1_lin_Norm-median_singlStart' & suffix = 'CR2082_HOLLOW_2Bands' 
  read_ldem,file,/ldem,/gauss1;,dir=dir

  
  
  
  ratio = sfbe/tfbe
  R     = total( abs(1.-ratio), 4 ) / float(nband)
  ZDA   = where(demc eq -999.)
  CNS   = where(demc ne -999.  AND R gt 0.25)

  Nesat = N_e
  Tmsat = Tm
  Rsat  = R


  superhigh=0.25 & p=where(demc ne -999. AND R ge superhigh) & if p(0) ne -1 then Rsat(p)=superhigh
  superlow=0.01  & p=where(demc ne -999. AND R le superlow)  & if p(0) ne -1 then Rsat(p)=superlow
  
  
  if ZDA(0) ne -1 then begin
     Nesat(ZDA)=0.
     Tmsat(ZDA)=0.
     Rsat (ZDA)=0.
  endif

  if CNS(0) ne -1 then begin
     Nesat(CNS)=0.
     Tmsat(CNS)=0.
     Rsat (CNS)=0.
  endif


  
  r0A=[1.025,1.065,1.105]


  xdisplay,map=Nesat,file='Ne_'+suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl=  4, titulo='Ne [10!U8!Ncm!U-3!N]' ,units=1.e8, minA=(r0A*0.+1.E-6), maxA=[2.2,1.6,1.0]  
  xdisplay,map=tmsat,file='Tm_'+suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl=  5, titulo='Tm [MK]'              ,units=1.e6, minA=(r0A*0.+1.E-6) ,maxA=[2.5,2.5,2.5]
  xdisplay,map=Rsat, file='R_' +suffix,nr=nr,nt=nth,rmin=1.0,rmax=1.3,r0A=r0A,win=0, clrtbl= 12, titulo='R'                               , minA=(r0A*0.+1.E-6) ,maxA=(r0A*0.+0.25)
  
  


  return
end
