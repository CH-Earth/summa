pro plot_paper2figure02

; define plotting parameters
window, 0, xs=1800, ys=700, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2.
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,2,0,1]

; define parameter set
iHRU = 1

; define constants
TFreeze = 273.16d

; define name of the NetCDF file
file_path = '/home/mclark/test_summa/summa/output/wrrPaperTestCases/figure02/'
file_name = file_path + 'vegImpactsWind_2006-2007_riparianAspenWindParamPerturb.nc'

; define the name of the validation data
valPath = '/home/mclark/test_summa/summa/testCases_data/validationData/'
valFile = valPath + 'ReynoldsCreek_eddyFlux.nc'

; define the path and name of the graphics file
gpath = '/home/mclark/test_summa/summa/verification/zFigures/'
gname = 'Clark_et_al__WRR2015b_figure02.png'

; *****
; * GET BASIC DATA FROM THE MODEL OUTPUT FILE...
; **********************************************

; open file
nc_file = ncdf_open(file_name, /nowrite)

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

 ; get the number of snow layers
 ivar_id = ncdf_varid(nc_file,'nSnow')
 ncdf_varget, nc_file, ivar_id, nSnow, offset=[iHRU,0], count=[1,ntime_mod]

 ; get the number of soil layers
 ivar_id = ncdf_varid(nc_file,'nSoil')
 ncdf_varget, nc_file, ivar_id, nSoil, offset=[iHRU,0], count=[1,ntime_mod]

 ; get the number of layers
 ivar_id = ncdf_varid(nc_file,'nLayers')
 ncdf_varget, nc_file, ivar_id, nLayers, offset=[iHRU,0], count=[1,ntime_mod]

 ; get the start index for midSoil
 ivar_id = ncdf_varid(nc_file,'midSoilStartIndex')
 ncdf_varget, nc_file, ivar_id, midSoilStartIndex, offset=[iHRU,0], count=[1,ntime_mod]

 ; get the start index for midToto
 ivar_id = ncdf_varid(nc_file,'midTotoStartIndex')
 ncdf_varget, nc_file, ivar_id, midTotoStartIndex, offset=[iHRU,0], count=[1,ntime_mod]

 ; get the start index for ifcToto
 ivar_id = ncdf_varid(nc_file,'ifcTotoStartIndex')
 ncdf_varget, nc_file, ivar_id, ifcTotoStartIndex, offset=[iHRU,0], count=[1,ntime_mod]

 ; get the start index for ifcToto
 ivar_id = ncdf_varid(nc_file,'ifcSnowStartIndex')
 ncdf_varget, nc_file, ivar_id, ifcSnowStartIndex, offset=[iHRU,0], count=[1,ntime_mod]

 ; define the date format
 dummy = label_date(date_format=['%D!C%M'])

 ;ibeg_mod = 0
 ;iend_mod = ntime_mod-1

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
 isubset  = where((iy_obs eq iy_mod[iend_mod]) and (im_obs eq im_mod[iend_mod]) and (id_obs eq id_mod[iend_mod]) and (ih_obs eq ih_mod[ibeg_mod]), nsubset)
 if(nsubset gt 0)then begin
  iend_obs = isubset[nsubset-1]
 endif else begin
  ;stop, 'no validation data in the simulation time period'
 endelse

; close NetCDF file
ncdf_close, nc_file


; *****
; * PLOT DIFFERENCE IN WIND SPEED...
; **********************************

print, '*** wind speed ***'

; define observations of interest
obsnames = ['wind', 'wind']
obs_stns = [1,       0]  ; 0 = Aspen understory; 1 = Aspen; 2 = Sagebrush

; define model variables of interest
modnames = ['windspd','scalarWindspdCanopyBottom']

; define y title
ytit = ['Above-canopy!Cwindspeed (m s!e-1!n)', 'Below-canopy!Cwindspeed (m s!e-1!n)']

; define legend
yleg1 = ['Observed windspeed above the Aspen canopy','Observed windspeed below the Aspen canopy']
yleg2 = ['Observed windspeed at the "exposed" site','Simulated windspeed below the Aspen canopy']

; loop through variables
for ivar=0,n_elements(modnames)-1 do begin

 ; get observations
 nc_file = ncdf_open(valFile, /nowrite)
  ivar_id = ncdf_varid(nc_file,obsnames[ivar])
  ncdf_varget, nc_file, ivar_id, obsdata, offset=[obs_stns[ivar],0], count=[1,ntime_obs]
 ncdf_close, nc_file

 ; get model output
 nc_file = ncdf_open(file_name, /nowrite)
  ivar_id = ncdf_varid(nc_file,modnames[ivar])
  ncdf_varget, nc_file, ivar_id, moddata, offset=[iHRU,0], count=[1,ntime_mod]
 ncdf_close, nc_file

 jtime=ibeg_obs
 obsmatch = dblarr(n_elements(moddata))
 ; match obs with the model
 for itime=0,n_elements(moddata)-1 do begin
  ; check we have matched correctly
  ;print, im_mod[itime], id_mod[itime], iy_mod[itime], ih_mod[itime], imi_mod[itime]
  ;print, im_obs[jtime], id_obs[jtime], iy_obs[jtime], ih_obs[jtime], imi_obs[jtime]
  ;print, im_obs[jtime+1], id_obs[jtime+1], iy_obs[jtime+1], ih_obs[jtime+1], imi_obs[jtime+1]
  if((iy_obs[jtime] ne iy_mod[itime]) or $
     (im_obs[jtime] ne im_mod[itime]) or $
     (id_obs[jtime] ne id_mod[itime]) or $
     (ih_obs[jtime] ne ih_mod[itime]) )then stop, 'dates do not match'
  if(ih_obs[jtime] ne ih_obs[jtime+1])then stop, 'different hour'
  ; average the wind speed
  wind = obsdata[jtime:jtime+1]
  if(min(wind) gt 0.d)then begin
   obsmatch[itime] = mean(wind)
  endif else begin
   obsmatch[itime] = -9999.d
  endelse
  ; increment index
  jtime=jtime+2
 endfor


 ; get hourly obs for the model
 ;for itime=0,n_elements(obsdata)-1, 2 do begin
 ; wind = obsdata[itime:itime+1]
 ; if(min(wind) gt 0.d)then begin
 ;  obsdata[itime:itime+1] = mean(wind)
 ; endif else begin
 ;  obsdata[itime:itime+1] = -9999.d
 ; endelse
 ;endfor

 if(ivar eq 0)then ymar = [3,1]
 if(ivar eq 1)then ymar = [4,0]

 ; make a base plot
 plot, djulian_mod[ibeg_mod:iend_mod], xrange=[djulian_mod[ibeg_mod],djulian_mod[iend_mod]], yrange=[0,20], xstyle=1, ystyle=1, $
  xtickformat=['label_date'], xticks=12, ytitle = ytit[ivar], xmargin=[70,-70], ymargin = ymar, title=ptitle, /nodata

 ; plot data
 ;oplot, djulian_obs[ibeg_obs:iend_obs], obsdata[ibeg_obs:iend_obs], min_value=(-9998.), color=80
 oplot, djulian_mod[ibeg_mod:iend_mod], obsmatch[ibeg_mod:iend_mod], min_value=(-9998.), color=0
 oplot, djulian_mod[ibeg_mod:iend_mod], moddata[ibeg_mod:iend_mod], color=210

 ; make a legend
 x1 = djulian_mod[ibeg_mod] + 0.48*(djulian_mod[iend_mod] - djulian_mod[ibeg_mod])
 x2 = djulian_mod[ibeg_mod] + 0.58*(djulian_mod[iend_mod] - djulian_mod[ibeg_mod])
 plots, [x1,x2], [18,18]
 plots, [x1,x2], [16,16], color=210
 xyouts, x2+3, 17.75, yleg1[ivar], charsize=1.5
 xyouts, x2+3, 15.75, yleg2[ivar], charsize=1.5

endfor 

; make a base plot
plot, indgen(5), xrange=[0,1], yrange=[0,5], xstyle=1, ystyle=1, $
 xtitle='Exceedance probability', ytitle='Below-canopy windspeed (m s!e-1!n)', $
 ymargin=[-13.5,1], xmargin=[-70,90], /nodata

; identify valid data
ivalid = where(obsmatch gt 0., nvalid)

; compute exceedance probability for the observations
windObs = obsmatch[ivalid]
sortObs = windObs(sort(windObs))
cumProb = dindgen(nvalid)/double(nvalid+1)
oplot, 1.d - cumProb, sortObs, color=0, thick=2

; define colors
icolor=[254,210,160,80,40]

; get different parameter sets
for iHRU=0,3 do begin

 ; get model output
 nc_file = ncdf_open(file_name, /nowrite)
  ivar_id = ncdf_varid(nc_file,'scalarWindspdCanopyBottom')
  ncdf_varget, nc_file, ivar_id, moddata, offset=[iHRU,0], count=[1,ntime_mod]
 ncdf_close, nc_file

 ; compute exceedance probability for the model
 windMod = moddata[ivalid]
 sortMod = windMod(sort(windMod))
 cumProb = dindgen(nvalid)/double(nvalid+1)
 oplot, 1.d - cumProb, sortMod, color=icolor[iHRU], thick=2

endfor ; looping through HRUs

xline = [0.4, 0.6]
y_top = 4.75
y_gap = 0.25
xcolr = [254,210,160,80,0]
xtext = ['Canopy wind parameter = 0.10',$
		 'Canopy wind parameter = 0.28',$
         'Canopy wind parameter = 0.50',$
         'Canopy wind parameter = 0.75',$
         'Below-canopy observations']
for ivar=0,4 do begin
 yline = y_top - float(ivar)*y_gap
 plots, xline, [yline,yline], color=xcolr[ivar], thick=2
 xyouts, xline[1]+0.01, yline-0.025, xtext[ivar], charsize=1.5
endfor

; make a figure
write_png, gpath+gname, tvrd(true=1)


stop

end
