; Wrapper to execute filter_ratio routine
pro DO_filter_ratio
 
; Root directory
; root_dir='/media/Data1/data1/'
  root_dir='/data1/'  

; FBEs directory:
  FBEdir = root_dir+'tomography/bindata/'

; TRFs directory:
  Qkldir = root_dir+'DATA/EUV_TRFs_and_PassBands/'

; FBEs files:
  datafiles= ['x_euviA_171_CR2082_hollow_NEW_Rmin1.00_Rmax1.30_IRmax1.25_30x90x180_BF4_L0.38',$
              'x_euviA_195_CR2082_hollow_NEW_Rmin1.00_Rmax1.30_IRmax1.25_30x90x180_BF4_L0.28',$
              'x_euviA_284_CR2082_hollow_NEW_Rmin1.00_Rmax1.30_IRmax1.25_30x90x180_BF4_L0.35']

; Tomographic grid parameters
  nr = 30 &   nt = 90 &   np = 180

; Bands to use:
  bandsindexes = [0,1]          ; [0,2] ; [1,2]

; String to make the file-name
  file_string='CR2082_EUVI171-195A_NewGrid'

; Directory to save the results  
  DEMdir = root_dir+'DATA/ldem_files/'

; TRFs files 
  file_ioneq = 'chianti.ioneq'  
  file_abund = 'sun_coronal_1992_feldman_ext.abund'
  filter     = 'S1a'
  qklfiles   = Qkldir+['Qkl_171_','Qkl_195_','Qkl_284_']+$
                file_ioneq+'_'+file_abund+$
                '.'+filter+'-Filter-photons-Abund-1e-3-ALL-withCONTINUUM_Ne1E08_C71.out'
; Temperature grid parameters
  Tmin =  0.90;0.75                   ; MK
  Tmax =  2.05;2.25                   ; MK
  DT   =  5.e-4                       ; MK 

; call to main routine
  filter_ratio,FBEdir,datafiles,qklfiles,DEMdir,file_string,nr,nt,np,Tmin,Tmax,DT,bandsindexes

  return
end


; MAIN-ROUTINE
pro filter_ratio,FBEdir,datafiles,qklfiles,DEMdir,file_string,nr,nt,np,Tmin,Tmax,DT,bandsindexes


; Select the 2 bands to do filter-ratio
  datafiles = datafiles (bandsindexes)  
  qklfiles  = qklfiles  (bandsindexes)
  nband     = n_elements(bandsindexes)
  if nband ne 2 then stop
  
  
; Load the FBEs in a 4-dimensional array 
  FBE_A  = fltarr(nr,nt,np,nband)
  tmp    = fltarr(nr,nt,np)
  for i=0,nband-1 do begin
     openr,1,FBEdir+datafiles(i)
     readu,1,tmp
     FBE_A(*,*,*,i)=tmp
     close,1
  endfor
; CORRECT FBE UNITS TO ACCOUNT FOR W* USED IN RECONSTRUCTION
  FBE_A = FBE_A / 0.0696
 
; Change Tmin and Tmax to K and calculete de number of temp. bins
  Tmin = Tmin*1.e6                   ; K
  Tmax = Tmax*1.e6                   ; K
  DT   = DT  *1.e6                   ; K 
  nTbins=LONG((Tmax-Tmin)/DT) ; Number of Temp bins
; Read the TRF files
  read_Qkl,qklfiles,logT,kernel,nband
; Interpol the TRF to a linear grid and create the temperature grid
  Interpol_Qkl,logT,kernel,nband,nTbins,Tmin,Tmax,Qkl,X 
; Temperature grid in K  
  T=x*Tmin
; TRF curves for the two selected bands  
  Qa=reform(Qkl(0,*))
  Qb=reform(Qkl(1,*))

; Initialize the density and temperatur arrays
  Ne_ratio = fltarr(nr,nt,np)-999.
  Te_ratio = fltarr(nr,nt,np)-999.

  for ir=0,nr-1 do begin
     for it=0,nt-1 do begin
        for ip=0,np-1 do begin
        ;  Filter-Ratio in a voxel   
           Ia = FBE_A(ir,it,ip,0)
           Ib = FBE_A(ir,it,ip,1)
           if Ia gt 0. and Ib gt 0. then begin
              if Ia/Ib ge min(Qa/Qb) and Ia/Ib le max(Qa/Qb) then begin 
                 f  = abs(Qa/Qb-Ia/Ib)
                 p  = (where(f eq min(f)))(0)
                 T0 = T(p)
                 N0 = sqrt(mean([Ia/Qa(p),Ib/Qb(p)]))
              endif else begin
                 T0 = -666.
                 N0 = -666.
              endelse
              Ne_ratio(ir,it,ip) = N0
              Te_ratio(ir,it,ip) = T0
           endif
        endfor
     endfor
  endfor

; Test that arrays are float (important to read with xread routine)  
; help,Ne_ratio,Te_ratio
; STOP
  
; Save the results  
  openw,1,DEMdir+'Ne_ratio_'+file_string+'.dat'
  writeu,1,Ne_ratio
  close,1
  openw,1,DEMdir+'Te_ratio_'+file_string+'.dat'
  writeu,1,Te_ratio
  close,1

  return
end


; ======================= SUB-ROUTINES ==================================

pro read_Qkl,qklfiles,logT,kernel,nband

; Initialize variables to read  
  x        =''
  N_e      =0.
  logem    =0.
  inst_fwhm=0.
  Ntemp    =0.
  w1       =0.
  w2       =0.
  wbin     =0.
  Nlambda  =0.
  psi      =0.
  psi_const=fltarr(nband)
  tt       =0.
  ii       =0.

; Read the TRFs from ASCII files  
  for ip=0,nband-1 do begin ; band loop
     openr,1,qklfiles(ip)
     for i=1,4 do readf,1,x
     readf,1,n_e,logem
     readf,1,x
     readf,1,Inst_fwhm, NTemp, w1, w2, wbin, Nlambda
     readf,1,x
     readf,1,psi
     psi_const(ip)=psi
     readf,1,x
     if ip eq 0 then begin
        kernel=fltarr(nband,ntemp) ; Initialize Q_k array
        logT  =fltarr(ntemp)       ; Initialize log(T_e) array
     endif
     for i=0,ntemp-1 do begin   ; Temperature loop
        readf,1,tt,ii
        if ip eq 0 then logT(i)=tt
        kernel(ip,i)=ii
     endfor                     ; Finish Temperature loop
     close,1
  endfor                        ; Finish band loop
  kernel=kernel*1.d0            ; convert to double presition
  return
end


pro Interpol_Qkl,logT,kernel,nband,nTbins,Tmin,Tmax,Qkl,X

; Select range of interest
  range=where( (logT ge alog10(Tmin)) AND (logT le alog10(Tmax)) )
  range=[min(range)-1,range,max(range)+1]

; Select the temperature range to be used:
; and correct Tmin and Tmax
  logaT  = [logT(    range)]
  skernel = kernel(*,range)
; Tmin   = 10.^min(logaT) 
; Tmax   = 10.^max(logaT)

; derive x = T/Tmin from log(T)   
  Xold = 10.^logaT/Tmin

; Decide NEW grid
; Make  X = T/Tmin
  DT_L = (Tmax-Tmin) / float(nTbins)
  T_L  = Tmin + DT_L /2. + DT_L *findgen(nTbins)
  X    = T_L /Tmin
     
; Interpolate logaT into the X-grid
  logaT = interpol(logaT,Xold,X)

; Make the Qkl array, interpolated into the X-grid
  Qkl=dblarr(nband,nTbins)
  for ib=0,nband-1 do Qkl(ib,*) = interpol(skernel(ib,*),Xold,X) ;,/spline)

;Check interpolation:
;goto,skip
    window,0
    plot,10.^logT/1.e6,kernel(0,*),xr=[0.5,5],/ylog,xtitle='T!de!n [MK]',title='TRFs'
    for ib=0,nband-1 do oplot,10.^logT/1.e6   ,kernel(ib,*)
    for ib=0,nband-1 do oplot,X*Tmin/1.e6,        Qkl(ib,*),psym=4,th=4
    window,1
    if nband ne 2 then STOP
    plot,X*Tmin/1.e6,Qkl(0,*)/Qkl(1,*),xtitle='T!de!n [MK]',title='Ratio Qa/Qb'
    STOP
skip:
return
end
