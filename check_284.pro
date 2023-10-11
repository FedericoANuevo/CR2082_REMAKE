

pro check_284

  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

  dir  ='/media/Data1/data1/DATA/ldem_files/'
  file ='LDEM._CR2082_euvi.A_Hollow_3Bands_gauss1_lin_Norm-median_singlStart'  
  read_ldem,file,/ldem,/gauss1,dir=dir
  ratio = sfbe/tfbe
  R     = total( abs(1.-ratio), 4 ) / float(nband)
  ZDA   = where(demc eq -999.)
  CNS   = where(demc ne -999.  AND R gt 0.25)
  Irmin = 1.02
  Irmax = 1.25
  fbe171 = reform(tfbe(*,*,*,0))
  fbe195 = reform(tfbe(*,*,*,1))
  fbe284 = reform(tfbe(*,*,*,2))

  grid_array3D,rad,lat,lon,radA,latA,lonA

  index284 = where ( fbe171 gt 0. and fbe195 gt 0. and fbe284 le 0. and radA gt Irmin and radA lt Irmax and abs(latA) gt 50.)
  index    = where (                                                    radA gt Irmin and radA lt Irmax and abs(latA) gt 50.)



  C_k = fltarr(nband)
  for k=0,nband-1 do begin
     p      = where  ( reform(tfbe(*,*,*,k)) gt 0.)
     C_k(k) = median ( (tfbe(*,*,*,k))(p)         )  
  endfor
  
 
  STOP
  return
end
