pro xdisplay_fbes

; Path to use x-tools in SolarTom_idl  
  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH
  
; Graph settings
  clrtbl=34

; Instrumental rmin and rmax
  Irmin =  1.02
  Irmax =  1.25

; Grid settings
  rmin = 1.0
  rmax = 1.3
  nr   = 30
  nt   = 90
  np   = 2*nt
; Radial grid 
  dr = (rmax-rmin)/nr
  r  = rmin + dr*findgen(Nr) +dr/2
  irad = where ( r ge Irmin and r le Irmax)


; Root directory
; root_dir='/media/Data1/data1/'
  root_dir='/data1/'  
; FBE directory:
  FBEdir = root_dir+'tomography/bindata/'
; FBES files:
  datafiles= ['x_euviA_171_CR2082_hollow_NEW_Rmin1.00_Rmax1.30_IRmax1.25_30x90x180_BF4_L0.38',$
              'x_euviA_195_CR2082_hollow_NEW_Rmin1.00_Rmax1.30_IRmax1.25_30x90x180_BF4_L0.28',$
              'x_euviA_284_CR2082_hollow_NEW_Rmin1.00_Rmax1.30_IRmax1.25_30x90x180_BF4_L0.35']

  datafiles_fess = datafiles+'_FESS'

    
; sample heights to plot
  r0A=[1.075];[1.035,1.105]

; Read the FBE data-cube
; Using callsolve_cg
  xread,file=datafiles(0),nr=nr,nt=nt,np=np,map=x171 
  xread,file=datafiles(1),nr=nr,nt=nt,np=np,map=x195
  xread,file=datafiles(2),nr=nr,nt=nt,np=np,map=x284
; Using callsolve_fess
  xread,file=datafiles_fess(0),nr=nr,nt=nt,np=np,map=x171f 
  xread,file=datafiles_fess(1),nr=nr,nt=nt,np=np,map=x195f
  xread,file=datafiles_fess(2),nr=nr,nt=nt,np=np,map=x284f
  
; Check ZDAs
  print,'CALLSOLVE_CG:'
  tmp = x171(irad,*,*) & p = where (tmp le 0.)
  if p(0) ne -1 then print,'ZDAs in 171 A:',float(n_elements(p))/n_elements(tmp) 
  tmp = x195(irad,*,*) & p = where (tmp le 0.)
  if p(0) ne -1 then print,'ZDAs in 195 A:',float(n_elements(p))/n_elements(tmp) 
  tmp = x284(irad,*,*) & p = where (tmp le 0.)
  if p(0) ne -1 then print,'ZDAs in 284 A:',float(n_elements(p))/n_elements(tmp) 
  print,'CALLSOLVE_FESS:'
  tmp = x171f(irad,*,*) & p = where (tmp EQ 0.)
  if p(0) ne -1 then print,'ZDAs in 171 A:',float(n_elements(p))/n_elements(tmp) 
  tmp = x195f(irad,*,*) & p = where (tmp EQ 0.)
  if p(0) ne -1 then print,'ZDAs in 195 A:',float(n_elements(p))/n_elements(tmp) 
  tmp = x284f(irad,*,*) & p = where (tmp EQ 0.)
  if p(0) ne -1 then print,'ZDAs in 284 A:',float(n_elements(p))/n_elements(tmp) 
  STOP


  xdisplay,map=x171,clrtbl=clrtbl,$
           file='x171_CR2082',nr=nr,nt=nt,rmin=rmin,rmax=rmax,$
           r0A=r0A,win=2,titulo='FBE 171 A',minS=minS,maxS=maxS,/add_bw

  xdisplay,map=x171f,clrtbl=clrtbl,$
           file='x171_CR2082_FESS',nr=nr,nt=nt,rmin=rmin,rmax=rmax,$
           r0A=r0A,win=2,titulo='FBE 171 A (FESS)',minA=minS,maxA=maxS,/add_bw
  

   xdisplay,map=x195,clrtbl=clrtbl,$
           file='x195_CR2082',nr=nr,nt=nt,rmin=rmin,rmax=rmax,$
           r0A=r0A,win=2,titulo='FBE 195 A',minS=minS,maxS=maxS,/add_bw

   xdisplay,map=x195f,clrtbl=clrtbl,$
            file='x195_CR2082_FESS',nr=nr,nt=nt,rmin=rmin,rmax=rmax,$
            r0A=r0A,win=2,titulo='FBE 195 A (FESS)',minA=minS,maxA=maxS,/add_bw


   xdisplay,map=x284,clrtbl=clrtbl,$
           file='x284_CR2082',nr=nr,nt=nt,rmin=rmin,rmax=rmax,$
           r0A=r0A,win=2,titulo='FBE 284 A',minS=minS,maxS=maxS,/add_bw
 
   xdisplay,map=x284f,clrtbl=clrtbl,$
           file='x284_CR2082_FESS',nr=nr,nt=nt,rmin=rmin,rmax=rmax,$
           r0A=r0A,win=2,titulo='FBE 284 A (FESS)',minA=minS,maxA=maxS,/add_bw
 

  return
end
