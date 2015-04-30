pro daily2hourly_daymet

; used to make hourly forcing data for all sub-basins above tollgate

; define the data path
data_path = '/home/mclark/summa/input/tollgate/'

; define basin filename
filenm_bas = data_path + 'tollgateBasin_daily.nc'

; define station filename
filenm_stn = data_path + 'RME_forcing.nc'

; define desired variables in the station data
cDesire = ['pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum']

; define the hourly time step
dt_hour = 3600.d

; define the index of the desired site
iStn = 0

; open files
nc_bas = ncdf_open(filenm_bas, /nowrite)
nc_stn = ncdf_open(filenm_stn, /nowrite)

 ; get the time vectors
 djulian_bas = call_function('exTime',nc_bas)
 djulian_stn = call_function('exTime',nc_stn)

 ; get the year/month/day of the station data
 ; NOTE small offset to ensure hour 0 is the first day
 caldat, djulian_stn-0.00001d, jm, jd, jyyy, jh, jmin

 ; get the basin vector
 ivar_id = ncdf_varid(nc_bas,'basin')
 ncdf_varget, nc_bas, ivar_id, idBasin

 ; get the basin elevation
 ivar_id = ncdf_varid(nc_bas,'basinElev')
 ncdf_varget, nc_bas, ivar_id, basinElev

 ; compute atmospheric pressure (Pa)
 presSpatial = 101325.d * ( (293.d - 0.0065d * basinElev) / 293.d )^5.256d

 ; loop through basins
 for iBas=0,n_elements(idBasin)-1 do begin

  ; define ASCII file for the basin data
  filenm_ASCII = data_path + 'hourlyAscii/basinForcing.' + strtrim(string(idBasin[iBas],format='(i4.4)'),2) + '.txt'
  print, 'filenm_ASCII = ', filenm_ASCII

  ; open up the ASCII file for writing
  openw, out_unit, filenm_ASCII, /get_lun

  ; loop through the basin time series
  for iTime=1369,n_elements(djulian_bas)-1 do begin

   ; get the time for the basin data
   caldat, djulian_bas[iTime], im, id, iyyy
   if(iyyy eq 2008 and im eq 10)then break

   ; identify overlap with the station data
   iMatch = where(iyyy eq jyyy and im eq jm and id eq jd, nMatch)
   if(nMatch gt 0)then begin

    ; check
    ;print, im, id, iyyy, nMatch
    if(nMatch ne 24)then stop, 'expect 24 matches'

    ; loop through the variables
    for ivar=0,n_elements(cDesire)-1 do begin

     ; read the  station data subset
     ivar_id = ncdf_varid(nc_stn,cDesire[ivar])
     ncdf_varget, nc_stn, ivar_id, stnData, offset=[iStn,iMatch[0]], count=[1,24]

     ; identify variable
     case cDesire[ivar] of

      ; *
      ; ** process precipitation...
      ; ***************************
      'pptrate': begin

       ; read the total daily precipitation
       ivar_id = ncdf_varid(nc_bas,'prcp')
       ncdf_varget, nc_bas, ivar_id, prcp, offset=[iBas,iTime], count=[1,1]

       ; get the daily total
       xTotal = total(stnData)
       ;print, 'prcp, xTotal = ', prcp, xTotal

       ; get the ratio of hourly precip to daily total
       if(xTotal*86400.d gt 0.01d)then begin  ; > 0.01 mm/day
        xRatio = reform(stnData)/xTotal
       endif else begin
        xRatio = replicate(1.d/24.d, 24)
       endelse

       ; get the basin sw radiation
       prcpBasin = xRatio*prcp/dt_hour
       ;print, 'prcpBasin = ', prcpBasin

      end  ; precip

      ; *
      ; ** process sw radiation...
      ; **************************
      'SWRadAtm': begin

       ; read the average sw radiation
       ivar_id = ncdf_varid(nc_bas,'srad')
       ncdf_varget, nc_bas, ivar_id, srad, offset=[iBas,iTime], count=[1,1]

       ; read the day length
       ivar_id = ncdf_varid(nc_bas,'dayl')
       ncdf_varget, nc_bas, ivar_id, dayl, offset=[iBas,iTime], count=[1,1]

       ; get the total radiation (J m-2)
       tRad = sRad*dayl

       ; get the ratio of hourly radiation to daily total
       xRatio = reform(stnData)/total(stnData)

       ; get the basin sw radiation
       sRadBasin = xRatio*tRad/dt_hour
       ;print, 'sRadBasin = ', sRadBasin, format='(a,1x,24(f7.0,1x))'

      end  ; sw radiation

      ; *
      ; ** process sw radiation...
      ; **************************
      'LWRadAtm': begin

       ; just do a direct copy for now (no DayMet data)
       lRadBasin = reform(stnData)
       ;print, 'lRadBasin = ', lRadBasin, format='(a,1x,24(f7.0,1x))'

      end  ; lw radiation

      ; *
      ; ** process air temperature...
      ; *****************************
      'airtemp': begin

       ; read in maximum temperature
       ivar_id = ncdf_varid(nc_bas,'tmax')
       ncdf_varget, nc_bas, ivar_id, tmax, offset=[iBas,iTime], count=[1,1]

       ; read in minimum temperature
       ivar_id = ncdf_varid(nc_bas,'tmin')
       ncdf_varget, nc_bas, ivar_id, tmin, offset=[iBas,iTime], count=[1,1]

       ; get the average temperature (deg C --> K)
       tavg = 0.5d*(tmax + tmin) + 273.16
       ;print, 'tmax, tmin, tavg = ', tmax, tmin, tavg

       ; get the temperature anomalies
       xDiff = reform(stnData) - mean(stnData)
       ;print, 'xDiff     = ', xDiff, format='(a,1x,24(f7.4,1x))'

       ; add the anomalies to the mean value
       tavgBasin = tavg + xDiff
       ;print, 'tavgBasin = ', tavgBasin, format='(a,1x,24(f7.2,1x))'

      end  ; average temperature

      ; *
      ; ** process wind speed...
      ; ************************
      'windspd': begin

       ; just do a direct copy for now (no DayMet data)
       windBasin = reform(stnData)
       ;print, 'windBasin = ', windBasin, format='(a,1x,24(f7.4,1x))'

      end  ; wind speed

      ; *
      ; ** process air pressure...
      ; *******************************
      'airpres': begin

       ; just use temporal constant (function of elevation)
       presBasin = replicate(presSpatial[ibas],24)
       ;print, 'presBasin = ', presBasin, format='(a,1x,24(f7.0,1x))'

      end  ; air pressure

      ; *
      ; ** process specific humidity...
      ; *******************************
      'spechum': begin

       ; read in vapor pressure (Pa)
       ivar_id = ncdf_varid(nc_bas,'vp')
       ncdf_varget, nc_bas, ivar_id, vp, offset=[iBas,iTime], count=[1,1]

       ; compute the saturated vapor pressure (Pa)
       svpFrz = 610.8           ; Saturation water vapour pressure at 273.16K (Pa)
       tavgC  = tavg - 273.16d  ; temperature (degrees C)
       svp = svpFrz * exp( (17.27d * tavgC) / (237.30d + tavgC) ) ; Saturated Vapour Press (Pa)

       ; compute the relative humidity
       ; NOTE: rh constant over the day
       rh = min([vp/svp, 1.d])*100.d

       ; compute the specific humidity for each hour
       shumBasin =dblarr(24)
       for iHour=0,23 do begin
        td = call_function('RLHUM2DEWPT',tavgBasin[iHour],rh) ; dewpoint temperature
        shumBasin[iHour] = call_function('DEWPT2SPHM',td,presBasin[iHour])
       endfor
       ;print, 'shumBasin = ', shumBasin, format='(a,1x,24(f7.4,1x))'

      end  ; specific humidity

      ; *
      ; ** check that we got everything...
      ; **********************************
      else: stop, 'unable to identify the desired variable'

     endcase  ; idenitifying the desired variable
    
    endfor  ; looping through the variables

    ; write data
    for iHour=0,23 do begin
     vData = [prcpBasin[iHour], sradBasin[iHour], lradBasin[iHour], tavgBasin[iHour], windBasin[iHour], presBasin[iHour], shumBasin[iHour]]
     ;print, iyyy, im, id, iHour, 0, 0.d, vData, format='(i4,1x,4(i2,1x),f6.1,1x,e14.4,1x,5(f10.3,1x),e12.3)'
     printf, out_unit, iyyy, im, id, iHour, 0, 0.d, vData, format='(i4,1x,4(i2,1x),f6.1,1x,e14.4,1x,5(f10.3,1x),e12.3)'
    endfor

   endif  ; if there is a data match

  endfor  ; looping through the basin time series

  ; free up the ASCII unit
  free_lun, out_unit
  ;stop, 'finished basin'

 endfor  ; looping through the basins

; close netcdf files
ncdf_close, nc_bas
ncdf_close, nc_stn

stop
end


function exTime, nc_file

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
 if(strtrim(tunit_words[0],2) eq 'days')    then aoff=1.d
 if(strtrim(tunit_words[0],2) eq 'seconds') then aoff=1.d/86400.d

 ; extract the time vector
 ncdf_varget, nc_file, ivar_id, atime
 return, bjulian + atime*aoff

end



FUNCTION RLHUM2DEWPT, T, RLHUM
; Compute Dewpoint temperature from Relative Humidity
; ---- This is done with respect to water ONLY ----
;
; All units are SI standard - i.e. Kelvin and pascals
; Based on Tetens' formula (1930)
; Units note :              Pa = N m-2 = kg m-1 s-2
SATVPFRZ=     610.8       ; Saturation water vapour pressure at 273.16K (Pa)
W_RATIO =       0.622     ; molecular weight ratio of water to dry air (-)
TFREEZE =     273.16      ; freezing point of water

VPSAT = SATVPFRZ * EXP( (17.27*(T-TFREEZE)) / (237.30 + (T-TFREEZE)) ) ; sat vapor press at grid cell (Pa)
TDCEL = 237.30 * ALOG( (VPSAT/SATVPFRZ)*(RLHUM/100.) ) / $              ; dewpoint temperature         (C)
        (17.27 - ALOG( (VPSAT/SATVPFRZ)*(RLHUM/100.) ) )
return, TDCEL + TFREEZE

end

FUNCTION DEWPT2SPHM, DEWPT, PRESS
; Compute specific humidity from dewpoint temp with respect to water
; ---- This is done with respect to water ONLY ----
;
; All units are SI standard - i.e. Kelvin and pascals
; Based on Tetens' formula (1930)
; VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
; Units note :              Pa = N m-2 = kg m-1 s-2
SATVPFRZ=     610.8       ; Saturation water vapour pressure at 273.16K (Pa)
W_RATIO =       0.622     ; molecular weight ratio of water to dry air (-)
TFREEZE =     273.16      ; freezing point of water

TDCEL = DEWPT-TFREEZE
VPAIR = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )        ; Vapour Press           (Pa)
return, (VPAIR * W_RATIO)/(PRESS - (1.-W_RATIO)*VPAIR)       ; Specific humidity (g/g)

END
