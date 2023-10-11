pro plot_histo

; Path to use de x-tools  
  !PATH = Expand_Path('+/data1/tomography/SolarTom_idl') + ':' + !PATH

; Grid parameters
  nr   = 30
  nt   = 90
  np   = 2*nt
  rmin = 1.0
  rmax = 1.3
  dr   = (rmax-rmin)/nr

; directory and file-names 
  dir         = '/data1/DATA/ldem_files/'
  filesuffix1 = 'LDEM._CR2082_euvi.A_Hollow_3Bands_gauss1_lin_Norm-median_singlStart' 
  filesuffix2 = 'CR2082_HOLLOW_compound1.dat'
  filesuffix3 = 'CR2082_HOLLOW_compound2.dat'  
  file1_Ne    = 'Ne_'+filesuffix1 
  file1_Te    = 'Te_'+filesuffix1 
  file2_Ne    = 'Ne_'+filesuffix2 
  file2_Te    = 'Te_'+filesuffix2 
  file3_Ne    = 'Ne_'+filesuffix3 
  file3_Te    = 'Te_'+filesuffix3 

; Read de data-cubes  
  xread,dir=dir,file=file1_Ne,nr=nr,nt=nt,np=np,map=Ne1
  xread,dir=dir,file=file2_Ne,nr=nr,nt=nt,np=np,map=Ne2
  xread,dir=dir,file=file3_Ne,nr=nr,nt=nt,np=np,map=Ne3
  xread,dir=dir,file=file1_Te,nr=nr,nt=nt,np=np,map=Te1
  xread,dir=dir,file=file2_Te,nr=nr,nt=nt,np=np,map=Te2
  xread,dir=dir,file=file3_Te,nr=nr,nt=nt,np=np,map=Te3


; Directory where the graph are recorded
  dir_fig   = '/data1/tomography/SolarTom_idl/Figures/'
; X-title of histograms
  xtitle_Ne = 'N!de!n [10!U8!Ncm!U-3!N]'
  xtitle_Te = 'T!dm!n [MK]'
; Labels of histograms  
  label1    = '3band'
  label2    = 'comp1'
  label3    = 'comp2'
; vector of heights to analyze
  r0A       = [1.025,1.105,1.245]   
 
; NORTH CH histograms ===============================================================
  title_base= 'North CH'
  lat_range = [+60.,+90.] 
  lon_range = [  0.,360.]
  minA_Ne   = [0.5,0.2,0.05]
  maxA_Ne   = [1.5,0.6,0.20]
  suffix_ne = 'histo_Ne_CR2082_HOLLOW_NCH'
  suffix_te = 'histo_Te_CR2082_HOLLOW_NCH'
  for i=0,n_elements(r0A)-1 DO BEGIN
     height        = r0A[i]
     height_string = strmid(string(height),6,5)
     title = title_base+' at '+height_string+' R!Ds!N'
     height_string = (STRJOIN(STRSPLIT(height_string, /EXTRACT,'.'), '_'))
     fig_name_ne  = suffix_ne+'_'+height_string+'_Rsun'
     fig_name_te  = suffix_te+'_'+height_string+'_Rsun'
     rad_range = [height-dr/2,height+dr/2]
   ; Ne histograms
     plot_histo_3maps,Ne1/1.E8,Ne2/1.E8,Ne3/1.e8,rmin,rmax,nr,nt,np,title,$
                      xtitle_Ne,label1,label2,label3,$
                      dir_fig=dir_fig,fig_name=fig_name_ne,$
                      rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,$
                      xrange=[minA_Ne(i),maxA_Ne(i)]
   ; Te histograms
     plot_histo_3maps,Te1/1.E6,Te2/1.E6,Te3/1.e6,rmin,rmax,nr,nt,np,title,$
                      xtitle_Te,label1,label2,label3,$
                      dir_fig=dir_fig,fig_name=fig_name_te,$
                      rad_range=rad_range,lat_range=lat_range,lon_range=lon_range,$
                      xrange=[0.5,1.5]     
  ENDFOR
; NORTH CH histograms ===============================================================

  return
end
