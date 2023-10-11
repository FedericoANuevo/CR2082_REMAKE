pro xdisplay_filter_ratio



; Path to use x-tools in SolarTom_idl  
  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH
  
; Instrumental rmin and rmax
  Irmin =  1.02
  Irmax =  1.25

; Grid settings
  rmin = 1.0
  rmax = 1.3
  nr   = 30
  nt   = 90
  np   = 2*nt


  dirDEM   = '/data1/DATA/ldem_files/' 
; file_Ne  = 'Ne_ratio_CR2082_EUVI171-195A.dat'
; file_Te  = 'Te_ratio_CR2082_EUVI171-195A.dat'
  file_Ne  = 'Ne_ratio_CR2082_EUVI171-195A_NewGrid.dat'
  file_Te  = 'Te_ratio_CR2082_EUVI171-195A_NewGrid.dat'

  suffix='_CR2082_NG'

; Read the data-cube
  xread,file=file_Ne,nr=nr,nt=nt,np=np,map=Ne_ratio,dir=dirDEM 
  xread,file=file_Te,nr=nr,nt=nt,np=np,map=Te_ratio,dir=dirDEM

  index = where(Te_ratio gt 0.) 
  print,'min Te [MK]:',min(Te_ratio(index))/1.e6
  print,'max Te [MK]:',max(Te_ratio(index))/1.e6
  STOP


  r0A=[1.025,1.065,1.105,1.155,1.205,1.255]
  
  xdisplay,map=Ne_ratio,file='Ne_ratio'+suffix,nr=nr,nt=nt,rmin=rmin,rmax=rmax,r0A=r0A,clrtbl= 25,titulo='Ne [10!U8!Ncm!U-3!N]' ,units=1.e8, minA=(r0A*0.+1.E-6),$
           maxA=[2.0,1.75,1.5,1.25,1.0,0.75],/add_bw
  xdisplay,map=Te_ratio,file='Te_ratio'+suffix,nr=nr,nt=nt,rmin=rmin,rmax=rmax,r0A=r0A,clrtbl= 25,titulo='Tm [MK]'              ,units=1.e6, minA=(r0A*0.+1.E-6),$
           maxA=3.0+r0A*0.,/add_bw
  



  return
end
