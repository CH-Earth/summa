pro plot_paper2figure08

; define plotting parameters
window, 0, xs=2000, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,1,0,1]

; define constants
TFreeze = 273.16d
sb = 5.6705d-8

; define name of the NetCDF file
file_path = '/home/mclark/test_summa/summa/output/wrrPaperTestCases/'

; define file prefix for the soil stress experiments
file_pref0 = ['figure08/vegImpactsTranspire_','figure09/basinRunoff_']

; define the file suffix for the root profile experiments
file_suff1 = ['_perturbRoots','_perturbRoots','_perturbRoots']

; define the file suffix for the lateral flow experiments
file_suff2 = ['_1dRichards','_lumpedTopmodel','_distributedTopmodel']

; define the name of the validation data
valFile = 'zObs/ReynoldsCreek_eddyFlux.nc'

; define HRUs
ixHRU1 = [1,2,3]
ixHRU2 = [1,1,6]

; loop through years
for iyear=2006,2006 do begin

 ; get water year
 cWaterYear = strtrim(iyear,2)+'-'+strtrim(iyear+1,2)

 ; define the name of the netcdf file
 filenm = file_path + file_pref0[0] + cWaterYear + file_suff1[0] + '.nc'

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
  ntime_mod = n_elements(djulian_mod)

  ; define the date format
  dummy = label_date(date_format=['%D-%M!C%H:%I'])

  ; define desired range
  ;ibeg_mod = 0
  ;iend_mod = ntime_mod-1

  ; 30 May - 13 Aug
  ibeg_mod = 5784
  iend_mod = 7607

  ;ibeg_mod = 6311
  ;iend_mod = ibeg_mod + 24

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
  caldat, djulian_mod, im_mod, id_mod, iy_mod, ih_mod, imi_mod, asec
  caldat, djulian_obs, im_obs, id_obs, iy_obs, ih_obs, imi_obs, asec

  ; get start index
  isubset  = where((iy_obs eq iy_mod[ibeg_mod]) and (im_obs eq im_mod[ibeg_mod]) and (id_obs eq id_mod[ibeg_mod]) and (ih_obs eq ih_mod[ibeg_mod]), nsubset)
  if(nsubset gt 0)then begin
   ibeg_obs = isubset[0]
  endif else begin
   stop, 'no validation data in the simulation time period'
  endelse

  ; get end index 
  isubset  = where((iy_obs eq iy_mod[iend_mod]) and (im_obs eq im_mod[iend_mod]) and (id_obs eq id_mod[iend_mod]) and (ih_obs eq ih_mod[iend_mod]), nsubset)
  if(nsubset gt 0)then begin
   iend_obs = isubset[nsubset-1]
  endif else begin
   ;stop, 'no validation data in the simulation time period'
  endelse

 ; close NetCDF file
 ncdf_close, nc_file

 ; *****
 ; * GET USEFUL STUFF...
 ; *********************

 ; get dates
 caldat, djulian_mod[ibeg_mod:iend_mod], im_mod, id_mod, iy_mod, ih_mod, imi_mod

 ; get diurnal array
 datMod = dblarr(25)
 datObs = dblarr(49)

 ; define colors
 icolor=[80,160,210,250]

 ; get the xticks
 xticks = [' ',[strtrim(indgen(7)*3+3,2)],' ']

 ; define xtitle
 ytitle=['Total evapotranspiration (mm h!e-1!n)',' ',' ']

 ; loop through plots
 for iPlot=0,1 do begin

  ; *****
  ; * PLOT STATION DATA...
  ; **********************

  ; define station
  istn = 1  ; 0 = Aspen understory; 1 = Aspen; 2 = Sagebrush

  ; define variables of interest
  varname = 'LE-wpl'

  ymax = (-0.6d)

  ; define y margin
  if(iplot eq 0)then xmar=[14,2]
  if(iplot eq 1)then xmar=[8,8]

  ; make a base plot
  plot, indgen(24)+1, xrange=[0,24], yrange=[0,ymax], xstyle=9, ystyle=1, $
   xtitle='Time of day', xticks=8, xtickname=xticks, ytitle=ytitle[iplot], xcharsize=1.5, ycharsize=1.5, $
   xticklen=(-0.02), xmargin=xmar, ymargin=[6,2], /nodata
  plots, [0,24], [ymax,ymax]
 
  ; get undesirable wind directions
  nc_file = ncdf_open(valFile, /nowrite)
   ivar_id = ncdf_varid(nc_file, 'WindFlag')
   ncdf_varget, nc_file, ivar_id, iWindflag, offset=[istn,0], count=[1,ntime_obs]
  ncdf_close, nc_file
  iWindflagSubset = reform(iWindflag[ibeg_obs:iend_obs])

  ; get observed latent heat flux
  nc_file = ncdf_open(valFile, /nowrite)
   ivar_id = ncdf_varid(nc_file,varname)
   ncdf_varget, nc_file, ivar_id, vardata, offset=[istn,0], count=[1,ntime_obs]
  ncdf_close, nc_file
  obsSubset = reform(vardata[ibeg_obs:iend_obs])

  ; get observation times
  obsTime = ceil(ih_obs[ibeg_obs:iend_obs]*100+imi_obs[ibeg_obs:iend_obs]*1.66666666d)

  ; print the valid data
  ;for itime=ibeg_obs,iend_obs do begin
  ; print, iy_obs[itime], im_obs[itime], id_obs[itime], ih_obs[itime], vardata[0,itime], iWindflag[0,itime], $
  ;  format='(i4,1x,3(i2,1x),f13.5,1x,i4)'
  ;endfor

  jtime=0
  ixObsValid = bytarr(n_elements(obsSubset))
  ixModValid = bytarr(n_elements(iy_mod))
  ; match obs with the model
  for itime=0,n_elements(iy_mod)-1 do begin
   ; define time in date array
   ktime = jtime+ibeg_obs
   ; check we have matched correctly
   ;print, im_mod[itime], id_mod[itime], iy_mod[itime], ih_mod[itime], imi_mod[itime]
   ;print, im_obs[ktime], id_obs[ktime], iy_obs[ktime], ih_obs[ktime], imi_obs[ktime]
   ;print, im_obs[ktime+1], id_obs[ktime+1], iy_obs[ktime+1], ih_obs[ktime+1], imi_obs[ktime+1]
   ;stop
   if((iy_obs[ktime] ne iy_mod[itime]) or $
      (im_obs[ktime] ne im_mod[itime]) or $
      (id_obs[ktime] ne id_mod[itime]) or $
      (ih_obs[ktime] ne ih_mod[itime]) )then stop, 'dates do not match'
   if(ih_obs[ktime] ne ih_obs[ktime+1])then stop, 'different hour'
   ; check if the data is desired
   if(iWindflagSubset[jtime] eq 0 and iWindflagSubset[jtime+1] eq 0 and obsSubset[jtime] gt -9998.d and obsSubset[jtime+1] gt -9998.d)then begin
    ixObsValid[jtime:jtime+1] = 0
    ixModValid[itime]         = 0
   endif else begin
    ixObsValid[jtime:jtime+1] = 1
    ixModValid[itime]         = 1
   endelse
   ; increment index
   jtime=jtime+2
  endfor

  ; get 30-min averages over the desired period
  for jhour=0,47 do begin
   ;imatch = where(obsTime eq jhour*50 and iWindflagSubset eq 0 and obsSubset gt -9998.d, nmatch)
   imatch = where(obsTime eq jhour*50 and ixObsValid eq 0, nmatch)
   nrgAverage    = total(obsSubset[imatch])/double(nmatch)
   datObs[jhour] = (nrgAverage/2260000.d)*3600.d  ; W m-2 --> mm/hr
   ;print, jhour, datObs[jhour]
   ;print, obsSubset[imatch]
   ;stop
  endfor
  datObs[48] = datObs[0]  ; zee wrap-around

  ; get time and check
  xtime = dindgen(48)/2.d + 0.25d
  ;print, transpose([[xtime],[datObs[1:48]]])

  ; plot 30 min observations
  oplot, xtime, datObs[1:48], color=0, psym=sym(1), symsize=4

  ; *****
  ; * PLOT MODEL SIMULATIONS...
  ; ***************************

  ; define parameters
  parnames = ['theta_mp','theta_sat','theta_res','critSoilWilting','critSoilTranspire','vGn_alpha','vGn_n','mpExp','k_soil','k_macropore','rootingDepth','rootDistExp']

  ; loop through experiments
  for ifile=0,2 do begin

   ; define file prefix
   if(iplot le 1)then begin
    file_pref = file_pref0[iplot]
   endif else begin
    if(ifile eq 0)then begin 
     file_pref = file_pref0[0]
    endif else begin
     file_pref = file_pref0[iplot]
    endelse
   endelse

   ; define files
   case iplot of
    0: filenm = file_path + file_pref0[iplot] + cWaterYear + file_suff1[ifile] + '.nc'
    1: filenm = file_path + file_pref0[iplot] + cWaterYear + file_suff2[ifile] + '.nc'
   endcase

   ; define HRUs
   case iplot of
    0: iHRU = ixHRU1[ifile]-1
    1: iHRU = ixHRU2[ifile]-1
   endcase

   ; print progress
   print, filenm, iHRU

   ; get desired parameters
   nc_file = ncdf_open(filenm, /nowrite)
    for ipar=0,n_elements(parnames)-1 do begin
     ivar_id = ncdf_varid(nc_file,parnames[ipar])
     ncdf_varget, nc_file, ivar_id, xData
     print, parnames[ipar], xData, format='(a20,1x,100(f20.10,1x))'
    endfor
   ncdf_close, nc_file

   ; get simulated latent heat flux
   nc_file = ncdf_open(filenm, /nowrite)
    ivar_id = ncdf_varid(nc_file,'scalarLatHeatTotal')
    ncdf_varget, nc_file, ivar_id, vData, offset=[iHRU,0], count=[1,ntime_mod]
   ncdf_close, nc_file
   total_ET = (reform(vData[ibeg_mod:iend_mod])/2260000.d)*3600.d  ; W m-2 --> mm/hr

   ; get hourly averages of the model output over the desired period
   for jhour=0,23 do begin
    ;imatch = where(ih_mod eq jhour and ixModValid eq 0, nmatch)
    imatch = where(ih_mod eq jhour, nmatch)
    datMod[jhour] = total(total_ET[imatch])/double(nmatch)
   endfor
   datMod[24] = datMod[0]  ; zee wrap-around

   ; smooth with a 3-hour moving average
   datSmooth = fltarr(25)
   datSmooth[0]  = datMod[0]
   datSmooth[24] = datMod[24]
   for jhour=1,23 do begin
    datSmooth[jhour] = mean(datMod[jhour-1:jhour+1])
   endfor
 
   ; plot total ET
   xtime = dindgen(24)+0.5d
   ;print, transpose([[xtime],[datMod[1:24]]])
   oplot, xtime, datSmooth[1:24], color=icolor[ifile], thick=4

  endfor  ; (looping through experiments)

  ; define legend
  if(iplot eq 0)then xtext = ['Root exp = 1','Root exp = 0.5','Root exp = 0.25']
  if(iplot eq 1)then xtext = ["Baseflow = 1D Richards'",'Baseflow = Topmodel (lumped)','Baseflow = Topmodel (distributed)']

  ; plot legend
  xline = [1,3]
  y_top = -0.566
  y_gap = -0.033
  for jvar=0,2 do begin
   yline = y_top - float(jvar)*y_gap
   plots, xline, [yline,yline], color=icolor[jvar], thick=4
   xyouts, xline[1]+0.25, yline+0.005, xtext[jvar], charsize=3
  endfor
  plots, mean(xline), yline-y_gap, psym=sym(1), symsize=4
  xyouts, xline[1]+0.25, yline-y_gap+0.005, 'Observations', charsize=3

 endfor  ; (looping through plots)

endfor  ; (looping through years)


; write figure
write_png, 'zFigures/Clark_et_al__WRR2015b_figure08.png', tvrd(true=1)



stop


end
