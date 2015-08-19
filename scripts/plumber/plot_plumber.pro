pro plot_plumber

; define plotting parameters
window, 0, xs=2000, ys=1400, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=[0,4,5,0,1]

; define the path to the plumber data
plumber_path = '/d1/mclark/PLUMBER_data/'

; define the path to the model output
model_path = plumber_path + 'model_output/'

; define the path to the site data
site_path =  plumber_path + 'site_data/flux/'

; define the benchmark path
benchmark_path = plumber_path + 'benchmark_data/'

; define the model names
model_names = ['CABLE.2.0',                  $
               'CABLE_2.0_SLI.vxh599_r553',  $
               'CHTESSEL',                   $
               'COLASSiB.2.0',               $
               'ISBA_SURFEX_3l.SURFEX7.3',   $
               'ISBA_SURFEX_dif.SURFEX7.3',  $
               'JULES.3.1',                  $
               'JULES3.1_altP',              $
               'Mosaic.1',                   $
               'NOAH.2.7.1',                 $
               'Noah.3.2',                   $
               'NOAH.3.3',                   $
               'ORCHIDEE.trunk_r1401',       $
               'SUMMA.1.0'                   ]

; define the site names
site_names = ['Amplero',     $
              'Blodgett',    $
              'Bugac',       $
              'ElSaler2',    $
              'ElSaler',     $
              'Espirra',     $
              'FortPeck',    $
              'Harvard',     $
              'Hesse',       $
              'Howard',      $
              'Howlandm',    $
              'Hyytiala',    $
              'Kruger',      $
              'Loobos',      $
              'Merbleue',    $
              'Mopane',      $
              'Palang',      $
              'Sylvania',    $
              'Tumba',       $
              'UniMich'      ]

; define the number of models and sites
nModels = n_elements(model_names)
nSites  = n_elements(site_names)

; loop through sites
for iSite=0,nSites-1 do begin

 ; *****
 ; * READ IN THE SITE OBS...
 ; *************************

 ; define filename
 site_filename = site_path + site_names[iSite] + 'Fluxnet.1.4_flux.nc'

 ; open file for reading
 nc_file = ncdf_open(site_filename, /nowrite)

  ; get time units
  ivar_id = ncdf_varid(nc_file,'time')
  ncdf_attget, nc_file, ivar_id, 'units', bunits
  cunits = string(bunits)

  ; extract the units "words"
  tunit_words = strsplit(string(cunits),' ',/extract)
  tunit_idate = fix(strsplit(tunit_words[2],'-',/extract))
  tunit_ihour = fix(strsplit(tunit_words[3],':',/extract))
  bjulian     = julday(tunit_idate[1],tunit_idate[2],tunit_idate[0],tunit_ihour[0],tunit_ihour[1],tunit_ihour[2])

  ; get the offset in days
  if(strtrim(tunit_words[0],2) eq 'seconds') then aoff=1.d/86400.d else stop, 'unknown time units'

  ; extract the time vector
  ncdf_varget, nc_file, ivar_id, atime
  djulian = bjulian + atime*aoff

  ; get the number of time elements
  ntime = n_elements(djulian)-1

  ; read in sensible heat
  ivar_id = ncdf_varid(nc_file,'Qh')
  ncdf_varget, nc_file, ivar_id, Qh
  Qh = reform(Qh)

  ; read in sensible heat quality control flag
  ivar_id = ncdf_varid(nc_file,'Qh_qc')
  ncdf_varget, nc_file, ivar_id, Qh_qc
  Qh_qc = floor(reform(Qh_qc)+0.5)

  ; read in latent heat
  ivar_id = ncdf_varid(nc_file,'Qle')
  ncdf_varget, nc_file, ivar_id, Qle
  Qle = reform(Qle)

  ; read in latent heat quality control flag
  ivar_id = ncdf_varid(nc_file,'Qle_qc')
  ncdf_varget, nc_file, ivar_id, Qle_qc
  Qle_qc = floor(reform(Qle_qc)+0.5)

 ; close netcdf file
 ncdf_close, nc_file

 ; *****
 ; * READ IN THE BENCHMARK DATA...
 ; *******************************

 ; define filename
 benchmark_filename = benchmark_path + site_names[iSite] + 'Fluxnet_1.4_PLUMBER_benchmarks.nc'

 ; open file for reading
 nc_file = ncdf_open(benchmark_filename, /nowrite)

  ; read in linear regression for sensible heat
  ivar_id = ncdf_varid(nc_file,'Qh_1lin')
  ncdf_varget, nc_file, ivar_id, Qh_1lin
  Qh_1lin = reform(Qh_1lin)

  ; read in linear regression for latent heat
  ivar_id = ncdf_varid(nc_file,'Qle_1lin')
  ncdf_varget, nc_file, ivar_id, Qle_1lin
  Qle_1lin = reform(Qle_1lin)

 ; close netcdf file
 ncdf_close, nc_file

 ; *****
 ; * MAKE A BASE PLOT...
 ; *********************

 ; get dates
 caldat, djulian, im, id, iy, ih, imi, asec

 ; get observation times
 obsTime = ceil(ih*100+imi*1.66666666d)

 ; get the xticks
 xticks = [' ',[strtrim(indgen(7)*3+3,2)],' ']

 ; define ymax
 ymin = -350
 ymax =  350

 ; make the base plot
 plot, indgen(24)+1, xrange=[0,24], yrange=[ymin,ymax], xstyle=9, ystyle=1, $
   xmargin = [15,2], xticks=8, xtickname=xticks, ytitle='heat flux', $
   xcharsize=1.5, ycharsize=1.5, xticklen=(-0.02), title=site_names[iSite], $
   /nodata
 plots, [0,24], [0,0]
 plots, [0,24], [ymax,ymax]

 ; *****
 ; * GET THE MEAN DIURNAL CYCLES...
 ; ********************************

 ; get time of the day
 xtime = dindgen(48)/2.d + 0.25d

 ; get arrays for the mean diurnal cycle
 datQh  = fltarr(49)
 datQle = fltarr(49)

 ; get arrays for the linear regression
 linQh  = fltarr(49)
 linQle = fltarr(49)

 ; compute the mean diurnal cycle
 for jhour=0,47 do begin
  ;imatch = where(obsTime eq jhour*50 and Qh_qc eq 0 and Qle_qc eq 0, nmatch)
  imatch = where(obsTime eq jhour*50, nmatch)
  ; observations
  datQh[jhour]  = total(Qh[imatch])/float(nmatch)
  datQle[jhour] = total(Qle[imatch])/float(nmatch)
  ; model simulations
  ; linear regression
  linQh[jhour]  = total(Qh_1lin[imatch])/float(nmatch)
  linQle[jhour] = total(Qle_1lin[imatch])/float(nmatch)
 endfor  ; (looping through hours)

 ; wrap-around
 ; (obs)
 datQh[48]  = datQh[0]   ; zee wrap-around
 datQle[48] = datQle[0]  ; zee wrap-around
 ; (regress)
 linQh[48]  = linQh[0]   ; zee wrap-around
 linQle[48] = linQle[0]  ; zee wrap-around

 ; *****
 ; * GET THE MODEL SIMULATIONS...
 ; ******************************

 ; loop through models
 for imodel=0,nmodels-1 do begin

  ; define the model filename
  model_filename = model_path + model_names[imodel] + '/' + model_names[imodel] + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'
  print, imodel, ' - ', model_filename

  ; open file for reading
  nc_file = ncdf_open(model_filename, /nowrite)

   ; read in the sensible heat flux
   ivar_id = ncdf_varid(nc_file,'Qh')
   ncdf_varget, nc_file, ivar_id, Qh_mod
   Qh_mod = reform(Qh_mod)

   ; read in the latent heat flux
   ivar_id = ncdf_varid(nc_file,'Qle')
   ncdf_varget, nc_file, ivar_id, Qle_mod
   Qle_mod = reform(Qle_mod)

  ; close netcdf file
  ncdf_close, nc_file 

  ; change sign for CHTESSEL and SUMMA
  if(model_names[imodel] eq 'CHTESSEL' or model_names[imodel] eq 'SUMMA.1.0')then begin
   Qh_mod  = -Qh_mod
   Qle_mod = -Qle_mod
  endif

  ; get arrays for the model simulations
  modQh  = fltarr(49)
  modQle = fltarr(49)

  ; compute the mean diurnal cycle
  for jhour=0,47 do begin
   imatch = where(obsTime eq jhour*50, nmatch)
   modQh[jhour]  = total(Qh_mod[imatch])/float(nmatch)
   modQle[jhour] = total(Qle_mod[imatch])/float(nmatch)
  endfor  ; (looping through hours)
  modQh[48]  = modQh[0]   ; zee wrap-around
  modQle[48] = modQle[0]  ; zee wrap-around

  ; plot model simulations
  if(model_names[imodel] eq 'SUMMA.1.0')then begin
   oplot, xtime, modQh,  color=250, thick=2
   oplot, xtime, -modQle, color=60, thick=2
  endif else begin
   oplot, xtime, modQh,  color=210, thick=1
   oplot, xtime, -modQle, color=90, thick=1
  endelse

 endfor  ; (looping through models)

 ; plot observations
 oplot, xtime, datQh,  color=250, thick=4
 oplot, xtime, -datQle, color=60, thick=4

 ; plot the linear regression
 ;oplot, xtime, linQh,  color=250, thick=2
 ;oplot, xtime, -linQle, color=60, thick=2

endfor  ; (looping through sites)

; write figure
write_png, 'figures/plumber_results.png', tvrd(true=1)





stop
end
