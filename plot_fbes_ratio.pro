pro plot_fbes_ratio

  common comunes,tm,wt,nband,demc,PHI,parametrizacion,Tmin,Tmax,nr,nth,np,rad,lat,lon,lambda,WTc
  common results_tomo,tfbe,sfbe,N_e

  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

  rmin = 1.0
  rmax = 1.3
  r0A  = [1.025,1.065,1.105]

  dir  ='/data1/DATA/ldem_files/' ; para Home-Office
 ;dir  ='/media/Data1/data1/DATA/ldem_files/'

  file_3bands ='LDEM._CR2082_euvi.A_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' 
  read_ldem,file_3bands,/ldem,/gauss1,dir=dir

  
  fbe_171 = tfbe(*,*,*,0)
  fbe_195 = tfbe(*,*,*,1)
  fbe_284 = tfbe(*,*,*,2)

  index1 = where (fbe_171 le 0. or fbe_195 le 0.)
  index2 = where (fbe_171 le 0. or fbe_284 le 0.)
  index3 = where (fbe_195 le 0. or fbe_284 le 0.)

  ratio1 = fbe_195/fbe_171
  ratio2 = fbe_284/fbe_171
  ratio3 = fbe_284/fbe_195

  ratio1(index1)=0.
  ratio2(index2)=0.
  ratio3(index3)=0.


  xdisplay,map=ratio1, file='ratio_195_over_171',nr=nr,nt=nth,rmin=rmin,rmax=rmax,r0A=r0A,win=0, clrtbl= 12, titulo='Ratio 195/171',minA=1.E-6+r0A*0.,maxA=1.5+r0A*0.  
  xdisplay,map=ratio2, file='ratio_284_over_171',nr=nr,nt=nth,rmin=rmin,rmax=rmax,r0A=r0A,win=0, clrtbl= 12, titulo='Ratio 284/171',minA=1.E-6+r0A*0.,maxA=1.5+r0A*0.  
  xdisplay,map=ratio3, file='ratio_284_over_195',nr=nr,nt=nth,rmin=rmin,rmax=rmax,r0A=r0A,win=0, clrtbl= 12, titulo='Ratio 284/195',minA=1.E-6+r0A*0.,maxA=1.5+r0A*0.  

  return
end
