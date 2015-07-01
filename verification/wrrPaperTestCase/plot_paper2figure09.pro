pro plot_paper2figure09

; define plotting parameters
window, 1, xs=1800, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3
!P.COLOR=0
erase, color=255
;!P.MULTI=[0,1,3]
!P.MULTI=[0,1,2]

; define the date format
dummy = label_date(date_format=['%D-%M!C%Y'])

; define the color "grey"
tvlct, r, g, b, /get
r[1] = 180
g[1] = 180
b[1] = 180
tvlct, r, g, b

; define constants
Tfreeze = 273.16
iden_ice= 917.

; define basin area (m2)
basarea = 389700.0d ; from Pierson et al. (WRR2001)

; define some conversion constants
secprday = 86400.d  ; number of seconds in a day
iden_wat = 1000.d   ; intrinsic density of water (kg m-3)

; define the path and name of the graphics file
gpath = '/home/mclark/test_summa/summa/verification/zFigures/'
gname = 'Clark_et_al__WRR2015b_figure09.png'

; define the path and name of the validation data
valPath = '/home/mclark/test_summa/summa/testCases_data/validationData/'
valFile = valPath+'ReynoldsCreek_valData.nc'

; define variable name and conversion factor
valName = 'Q'
valMult = iden_wat*secprday/basarea

; define file path for the distributed simulations
file_path = '/home/mclark/test_summa/summa/output/wrrPaperTestCases/figure09/'

; define suffixes
cSuffix = ['_1dRichards','_distributedTopmodel']

; define prefix
fPrefix = 'basinRunoff_'

; define julian days
time0 = julday(10,1,2002)
time1 = julday(10,1,2008)

; define titles
ptit = ["Baseflow = 1D Richards'",'Baseflow = saturated sub-surface flow (distributed)']

; define colors
icolor=[80,250]

; define HRU for 1-d simulations
iHRU=0

; loop through different prefixes
for ifile=0,1 do begin

 ; define suffix and prefix
 iSuffix=ifile

 ; define y margin
 if(ifile eq 0)then ymar=[ 1,5]
 if(ifile eq 1)then ymar=[ 5,1]

 ; make a base plot for some time series
 if(ifile eq 0)then begin
  plot, indgen(5), xrange=[time0,time1], yrange=[0,35], xstyle=1, ystyle=1, xticklen=(-0.02),$
   xtickname=[' ',' ',' ',' ',' ',' ',' '], xticks=6, ytitle = 'Runoff!C(mm/day)', $
   xmargin=[12,5], ymargin=ymar, /nodata
 endif else begin
  plot, indgen(5), xrange=[time0,time1], yrange=[0,35], xstyle=1, ystyle=1, xticklen=(-0.02),$
   xtickformat=['label_date'], xticks=6, ytitle = 'Runoff!C(mm/day)', $
   xmargin=[12,5], ymargin=ymar, /nodata
 endelse

 ; plot the experiment
 xyouts, time0+50.d, 30, ptit[ifile], charsize=2.5

 ; plot the legend
 x1 = time0+100.d
 x2 = x1+200.d
 oplot, [x1,x2], [27.0,27.0], color=250
 xyouts, x2+5.d, 26.d, 'Simulations', charsize=2.

 ; define the observations
 y1 = 22.0d
 y2 = 24.5d
 polyfill, [x1,x2,x2,x1], [y1,y1,y2,y2], color=1
 xyouts, x2+5., y1, 'Observations', charsize=2.

 ; junk loop (not used)
 for iJunk=0,0 do begin

  ; loop through the years
  for iyear=2002,2007 do begin

   ; declare model output file
   filenm = file_path + fPrefix+strtrim(iyear,2)+'-'+strtrim(iyear+1,2)+cSuffix[ifile]+'.nc'
   print, filenm

   ; *****
   ; * GET THE VALIDATION DATA...
   ; ****************************

   ; open file
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

    ; get data subset
    caldat, djulian_obs, im, id, iy, ih, imi, asec
    isubset = where( (iy eq iyear and im ge 10) or (iy eq iyear+1 and im lt 10), nsubset)
    ibeg_obs = isubset[0]
    iend_obs = isubset[nsubset-1]

    ; get the desired variable
    ivar_id = ncdf_varid(nc_file,valName)
    ncdf_varget, nc_file, ivar_id, valData

    ; compute cumulative runoff
    ;cRunoff = dblarr(ntime)
    ;cRunoff[0] = valdata[0]*iden_wat*3600.d/basarea
    ;for itime=1,ntime-1 do begin
    ; cRunoff[itime] = cRunoff[itime-1] + valdata[itime+i_beg]*iden_wat*3600.d/basarea
    ;endfor
 
    ; plot runoff
    ;oplot, djulian[i_beg:i_end], cRunoff[0:ntime-1], color=80
 
   ; close the NetCDF file
   ncdf_close, nc_file


   ; *****
   ; * PLOT THE MODEL SIMULATIONS...
   ; *******************************

   ; open file
   nc_file = ncdf_open(filenm, /nowrite)

    ; define some parameters to check
    parlist = ['frozenPrecipMultip','k_macropore','theta_mp','theta_res','theta_sat','vGn_n','fieldCapacity','f_impede','k_soil','kAnisotropic','zScale_TOPMODEL','qSurfScale']

    ; get some parameters
    for ipar=0,n_elements(parlist)-1 do begin
     ; get parameter
     ivar_id = ncdf_varid(nc_file,parlist[ipar])
     ncdf_varget, nc_file, ivar_id, pardata
     ; print
     if(iyear eq 2002)then print, strtrim(parlist[ipar],2)+': ', pardata, format='(a,100(f12.8,1x))'
    endfor

    ; get the HRU area
    ivar_id = ncdf_varid(nc_file,'HRUarea')
    ncdf_varget, nc_file, ivar_id, HRUarea

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

    ; define the date format
    dummy = label_date(date_format=['%D-%M!C%Y'])

    ; get the number of time elements
    ntime = n_elements(djulian)-1

    ; get start and end indices
    i_beg = 0
    i_end = ntime-1

    ; get the name of the groundwater decision
    NCDF_ATTGET, nc_file, 'groundwatr', btext, /global

    ; get the runoff
    if(string(btext) eq 'noXplict')then begin
     ivar_id = ncdf_varid(nc_file,'averageSurfaceRunoff')
     ncdf_varget, nc_file, ivar_id, averageSurfaceRunoff, offset=[iHRU,0], count=[1,ntime]
     ivar_id = ncdf_varid(nc_file,'averageSoilDrainage')
     ncdf_varget, nc_file, ivar_id, averageSoilDrainage, offset=[iHRU,0], count=[1,ntime]
     averageInstantRunoff = reform(averageSoilDrainage) + reform(averageSurfaceRunoff)
    endif else begin
     ivar_id = ncdf_varid(nc_file,'averageInstantRunoff')
     ncdf_varget, nc_file, ivar_id, averageInstantRunoff, offset=[0], count=[ntime]
    endelse
 
    ; get total routed runoff
    ivar_id = ncdf_varid(nc_file,'averageRoutedRunoff')
    ncdf_varget, nc_file, ivar_id, averageRoutedRunoff, offset=[0], count=[ntime]

   ; close NetCDF file
   ncdf_close, nc_file

   ; validation data
   ;if(iSuffix eq 0)then begin
    for itime=ibeg_obs,iend_obs do begin
     oplot, [djulian_obs[itime],djulian_obs[itime]], [0.,valdata[itime]*valMult], color=1
    endfor
    ;oplot, djulian_obs[ibeg_obs:iend_obs], valdata[ibeg_obs:iend_obs]*valMult, color=80
   ;endif

   ; plot runoff
   if(string(btext) eq 'noXplict')then begin
    oplot, djulian[i_beg:i_end], averageInstantRunoff[i_beg:i_end]*secprday*iden_wat, color=250
   endif else begin
    oplot, djulian[i_beg:i_end], averageRoutedRunoff[i_beg:i_end]*secprday*iden_wat, color=250
   endelse

   ;stop, 'checking a given year'

   ; save data
   if(iyear eq 2002)then begin
    tData = djulian_obs[ibeg_obs:iend_obs]
    xData = valdata[ibeg_obs:iend_obs]*valMult
    yData = averageRoutedRunoff[i_beg:i_end]*secprday*iden_wat
   endif else begin
    tData = [tData,djulian_obs[ibeg_obs:iend_obs]]
    xData = [xData,valdata[ibeg_obs:iend_obs]*valMult]
    yData = [yData,averageRoutedRunoff[i_beg:i_end]*secprday*iden_wat]
   endelse

  endfor  ; loop through years
  ;stop

  ; compute NS
  xMean = mean(xData)
  yMean = mean(yData)
  xNash = 1.d - total( (xData - yData)^2. ) / total( (xData - xMean)^2. )
  ;xyouts, djulian[i_beg] -800.d, 30.d, 'NS = '+strtrim(string(xNash,format='(f9.3)'),2)
  print, 'xNash, xMean, yMean = ', xNash, xMean, yMean


 endfor ; loop through parameters

endfor ; loop through experiments

 ; make a figure
write_png, gpath+gname, tvrd(true=1)

stop
end
