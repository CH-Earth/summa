pro plot_paper2figure01

; define plotting parameters
window, 0, xs=1200, ys=800, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3.
!P.COLOR=0
erase, color=255
!P.MULTI=[0,3,2,0,1]

; define parameter set
ipar = 0

; define constants
TFreeze = 273.16d
sb = 5.6705d-8

; define name of the NetCDF file
file_path = '/home/mclark/test_summa/summa/output/wrrPaperTestCases/figure01/'
file_pref = 'vegImpactsRad_'
file_suff = ['_riparianAspenBeersLaw','_riparianAspenNLscatter','_riparianAspenUEB2stream','_riparianAspenCLM2stream']

; define the name of the validation data
valFile = 'zObs/ReynoldsCreek_eddyFlux.nc'

; loop through years
for iyear=2005,2007 do begin

 ; get water year
 cWaterYear = strtrim(iyear,2)+'-'+strtrim(iyear+1,2)

 ; define the name of the netcdf file
 filenm = file_path + file_pref + cWaterYear + file_suff[0] + '.nc'

 ; *****
 ; * GET BASIC DATA FROM THE MODEL OUTPUT FILE...
 ; **********************************************

 ; open file
 nc_file = ncdf_open(filenm, /nowrite)

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
  djulian_mod = bjulian + atime*aoff

  ; get the number of time elements
  ntime_mod = n_elements(djulian_mod)-1

  ; define the date format
  dummy = label_date(date_format=['%D-%M!C%H:%I'])

  ; define desired range
  ibeg_mod = 0
  iend_mod = ntime_mod-1

 ; close NetCDF file
 ncdf_close, nc_file

 ; *****
 ; * GET BASIC DATA FROM THE VALIDATION DATA...
 ; ********************************************

 ; open file for reading
 nc_file = ncdf_open(valFile, /nowrite)

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
  djulian_obs = bjulian + atime*aoff

  ; get number of time elements
  ntime_obs = n_elements(djulian_obs)

  ; get data subset
  caldat, djulian_mod, im_mod, id_mod, iy_mod, ih, imi, asec
  caldat, djulian_obs, im_obs, id_obs, iy_obs, ih, imi, asec

  ; get start index
  isubset  = where((iy_obs eq iy_mod[ibeg_mod]) and (im_obs eq im_mod[ibeg_mod]) and (id_obs eq id_mod[ibeg_mod]), nsubset)
  if(nsubset gt 0)then begin
   ibeg_obs = isubset[0]
  endif else begin
   stop, 'no validation data in the simulation time period'
  endelse

  ; get end index 
  isubset  = where((iy_obs eq iy_mod[iend_mod]) and (im_obs eq im_mod[iend_mod]) and (id_obs eq id_mod[iend_mod]), nsubset)
  if(nsubset gt 0)then begin
   iend_obs = isubset[nsubset-1]
  endif else begin
   ;stop, 'no validation data in the simulation time period'
  endelse

 ; close NetCDF file
 ncdf_close, nc_file

 ; *****
 ; * PLOT DIFFERENCE IN CANOPY SHORTWAVE...
 ; ****************************************

 ; define the number of valid days in each year
 cdays = ['209','236','220']

 ; get the xticks
 xticks = [' ',[strtrim(indgen(7)*3+3,2)],' ']

 ; define the plot title
 ptitle = cWaterYear + ' (valid days = ' + cdays[iyear-2005] +')'

 ; define ytitle
 if(cWaterYear eq '2005-2006')then ytitle='Below canopy!Cshortwave radiation (W m!e-2!n)'
 if(cWaterYear ne '2005-2006')then ytitle=' '

 ; define xmargin
 if(cWaterYear eq '2005-2006')then xmar=[14,-2]
 if(cWaterYear eq '2006-2007')then xmar=[10,2]
 if(cWaterYear eq '2007-2008')then xmar=[ 6,6]

 ; make a base plot
 plot, indgen(24)+1, xrange=[0,24], yrange=[0,1000], xstyle=9, ystyle=1, $
  xtitle=' ', xticks=8, xtickname=xticks, ytitle=ytitle, xcharsize=1.5, ycharsize=1.5, $
  xticklen=(-0.02), xmargin=xmar, ymargin = [4,4], title=pTitle, /nodata
 plots, [0,24], [1000,1000]

 ; plot legend
 if(cWaterYear eq '2005-2006')then begin
  xline = [1,5]
  y_top = 950
  y_gap = 50
  xcolr = [0,80,160,210,250]
  xtext = ["Above canopy forcing","Beer's law","NL 1999","MT 2012","Dick. 1983"]
  for ivar=0,4 do begin
   yline = y_top - float(ivar)*y_gap
   plots, xline, [yline,yline], color=xcolr[ivar], thick=2
   xyouts, xline[1]+0.25, yline-10, xtext[ivar], charsize=1.25
  endfor
  plots, mean(xline), yline-y_gap, psym=sym(1), symsize=2
  xyouts, xline[1]+0.25, yline-y_gap-10, 'Obs', charsize=1.5
 endif

 ; define observations of interest
 obsname = 'rsd'
 obs_stn =     0  ; 0 = Aspen understory; 1 = Aspen; 2 = Sagebrush

 ; get observations
 nc_file = ncdf_open(valFile, /nowrite)
  ivar_id = ncdf_varid(nc_file,obsname)
  ncdf_varget, nc_file, ivar_id, obsdata, offset=[obs_stn,0], count=[1,ntime_obs]
 ncdf_close, nc_file

 ; get sw down
 nc_file = ncdf_open(filenm, /nowrite)
  ivar_id = ncdf_varid(nc_file,'SWRadAtm')
  ncdf_varget, nc_file, ivar_id, SWRadAtm, offset=[ipar,0], count=[1,ntime_mod]
  SWRadAtm = reform(SWRadAtm)
 ncdf_close, nc_file

 ; get dates
 caldat, djulian_obs[ibeg_obs:iend_obs], im_obs, id_obs, iy_obs, ih_obs, imi_obs
 caldat, djulian_mod[ibeg_mod:iend_mod], im_mod, id_mod, iy_mod, ih_mod, imi_mod

 ; get subset of observations
 obsSubset = obsdata[ibeg_obs:iend_obs]

 ; get subset of forcing
 frcSubset = SWRadAtm[ibeg_mod:iend_mod]

 ; get number of time elements
 ntime_obs = n_elements(obsSubset)
 ntime_mod = n_elements(frcSubset)

 ; define arrays for missing data
 imissing_obs = bytarr(ntime_obs)
 imissing_mod = bytarr(ntime_mod)

 ; get arrays for model and observations
 datObs = dblarr(49)
 datTop = dblarr(25)
 datMod = dblarr(25)

 ; get observation times
 obsTime = ceil(ih_obs*100+imi_obs*1.66666666d)

 ; identify days with missing data
 for itime=0,ntime_obs-1 do begin
  if(obsSubset[itime] lt -100.d)then begin
   ; find misisng day in observations
   imatch = where(im_obs eq im_obs[itime] and id_obs eq id_obs[itime], nmatch)
   if(nmatch ne 48)then stop, 'expect to find 48 data points in observations'
   imissing_obs[imatch] = 1
   ; find misisng day in model output
   imatch = where(im_mod eq im_obs[itime] and id_mod eq id_obs[itime], nmatch)
   if(nmatch lt 23)then stop, 'expect to find 24 data points in model output'
   imissing_mod[imatch] = 1
  endif  ; if missing data
 endfor

 ; get 30-min averages over the desired period
 for jhour=0,47 do begin
  imatch = where(obsTime eq jhour*50 and imissing_obs eq 0, nmatch)
  print, nmatch
  if(nmatch gt 50)then begin
   datObs[jhour] = total(obsSubset[imatch])/double(nmatch)
  endif else begin
   datObs[jhour] = 0.d
  endelse
 endfor
 datObs[48] = datObs[0]  ; zee wrap-around

 ; plot 30 min observations
 oplot, dindgen(48)/2.d + 0.25, datObs[1:48], color=0, psym=sym(1), symsize=2

 ; get hourly averages of the forcing over the desired period
 for jhour=0,23 do begin
  imatch = where(ih_mod eq jhour and imissing_mod eq 0, nmatch)
  print, nmatch
  if(nmatch gt 50)then begin
   datTop[jhour] = total(frcSubset[imatch])/double(nmatch)
  endif else begin
   datTop[jhour] = 0.d
  endelse
 endfor
 datTop[24] = datTop[0]  ; zee wrap-around

 ; plot forcing
 oplot, dindgen(24)+0.5d, datTop[1:24], color=0, thick=2

 ; define colors
 icolor=[80,160,210,250]

 ; loop through the different experiments
 for isuff=0,n_elements(file_suff)-1 do begin

  ; define the name of the netcdf file
  filenm = file_path + file_pref + cWaterYear + file_suff[isuff] + '.nc'

  ; get below canopy solar
  nc_file = ncdf_open(filenm, /nowrite)
   ivar_id = ncdf_varid(nc_file,'scalarBelowCanopySolar')
   ncdf_varget, nc_file, ivar_id, scalarBelowCanopySolar, offset=[ipar,0], count=[1,ntime_mod]
   scalarBelowCanopySolar = reform(scalarBelowCanopySolar)
  ncdf_close, nc_file

  ; get subset of below-canopy solar
  modSubset = scalarBelowCanopySolar[ibeg_mod:iend_mod]

  ; get hourly averages of the forcing over the desired period
  for jhour=0,23 do begin
   imatch = where(ih_mod eq jhour and imissing_mod eq 0, nmatch)
   if(nmatch gt 50)then begin
    datMod[jhour] = total(modSubset[imatch])/double(nmatch)
   endif else begin
    datMod[jhour] = 0.d
   endelse
  endfor
  datMod[24] = datMod[0]  ; zee wrap-around

  ; plot below-canopy solar
  oplot, dindgen(24)+0.5d, datMod[1:24], color=icolor[isuff], thick=2

 endfor  ; (looping through files)

 ; *****
 ; * PLOT RESULTS FROM PARAMETER PERTURBATION EXPERIMNENTS...
 ; **********************************************************

 ; make a base plot
 plot, indgen(24)+1, xrange=[0,24], yrange=[0,1000], xstyle=9, ystyle=1, $
  xtitle='Time of day', xticks=8, xtickname=xticks, ytitle=ytitle, xcharsize=1.5, ycharsize=1.5, $
  xticklen=(-0.02), xmargin=xmar, ymargin = [8,0], /nodata
 plots, [0,24], [1000,1000]

 ; define colors
 icolor=[40,80,160,210,250]

 ; plot 30 min observations
 oplot, dindgen(48)/2.d + 0.25, datObs[1:48], color=0, psym=sym(1), symsize=2

 ; plot forcing
 oplot, dindgen(24)+0.5d, datTop[1:24], color=0, thick=2

 ; define the name of the netcdf file
 filenm = file_path + file_pref + cWaterYear + '_riparianAspenVegParamPerturb.nc'

 ; loop through HRUs
 for iHRU=0,4 do begin

  ; get below canopy solar
  nc_file = ncdf_open(filenm, /nowrite)
   ivar_id = ncdf_varid(nc_file,'scalarBelowCanopySolar')
   ncdf_varget, nc_file, ivar_id, scalarBelowCanopySolar, offset=[iHRU,0], count=[1,ntime_mod]
   scalarBelowCanopySolar = reform(scalarBelowCanopySolar)
  ncdf_close, nc_file

  ; get subset of below-canopy solar
  modSubset = scalarBelowCanopySolar[ibeg_mod:iend_mod]

  ; get hourly averages of the forcing over the desired period
  for jhour=0,23 do begin
   imatch = where(ih_mod eq jhour and imissing_mod eq 0, nmatch)
   if(nmatch gt 50)then begin
    datMod[jhour] = total(modSubset[imatch])/double(nmatch)
   endif else begin
    datMod[jhour] = 0.d
   endelse
  endfor
  datMod[24] = datMod[0]  ; zee wrap-around

  ; plot below-canopy solar
  oplot, dindgen(24)+0.5d, datMod[1:24], color=icolor[iHRU], thick=2

 endfor  ; (looping through HRUs)

 ; plot legend
 if(cWaterYear eq '2005-2006')then begin
  xline = [1,5]
  y_top = 950
  y_gap = 50
  xcolr = [0,40,80,160,210,250]
  xtext = ["Above canopy forcing","LAI = 1.00","LAI = 1.25","LAI = 1.50","LAI = 1.75","LAI = 2.00"]
  for ivar=0,5 do begin
   yline = y_top - float(ivar)*y_gap
   plots, xline, [yline,yline], color=xcolr[ivar], thick=2
   xyouts, xline[1]+0.25, yline-10, xtext[ivar], charsize=1.25
  endfor
  plots, mean(xline), yline-y_gap, psym=sym(1), symsize=2
  xyouts, xline[1]+0.25, yline-y_gap-10, 'Obs', charsize=1.5
 endif

endfor  ; (looping through years)

; write figure
write_png, 'zFigures/Clark_et_al__WRR2015b_figure01.png', tvrd(true=1)



stop


end
