pro stn2grid_tollgateData

; *****
; (0) PRELIMINARIES...
; ********************

; define station data path
stn_path = '/home/mclark/summa/input/tollgate/stationData/'

; define station data file
stn_file = stn_path + 'tollgate_forcing.nc'

; define gridded data path
grid_path = '/home/mclark/summa/input/tollgate/'

; define gridded data file
grid_file = grid_path + 'stn2grid_tollgate.nc'

; *****
; (1) READ IN THE DAYMET METADATA...
; **********************************

; open netcdf file for reading
file_id = ncdf_open(grid_file, /nowrite)

 ; read the latitude
 ivarid = ncdf_varid(file_id,'lat')
 ncdf_varget, file_id, ivarid, grid_lat

 ; read the longitude
 ivarid = ncdf_varid(file_id,'lon')
 ncdf_varget, file_id, ivarid, grid_lon

 ; read the elevation
 ivarid = ncdf_varid(file_id,'elev')
 ncdf_varget, file_id, ivarid, grid_elev

 ; read the mask
 ivarid = ncdf_varid(file_id,'mask')
 ncdf_varget, file_id, ivarid, grid_mask

 ; get time units
 ivar_id = ncdf_varid(file_id,'time')
 ncdf_attget, file_id, ivar_id, 'units', bunits
 cunits = string(bunits)

 ; get the base julian day
 tunit_words  = strsplit(string(cunits),' ',/extract)
 tunit_idate  = fix(strsplit(tunit_words[2],'-',/extract))
 tunit_ihour  = fix(strsplit(tunit_words[3],':',/extract))
 bjulian_grid = julday(tunit_idate[1],tunit_idate[2],tunit_idate[0],tunit_ihour[0],tunit_ihour[1],tunit_ihour[2])

; close the file
ncdf_close, file_id

; *****
; (2) READ IN THE STATION METADATA...
; ***********************************

; open netcdf file for reading
file_id = ncdf_open(stn_file, /nowrite)

 ; read the latitude
 ivarid = ncdf_varid(file_id,'latitude_wgs84')
 ncdf_varget, file_id, ivarid, stn_lat

 ; read the longitude
 ivarid = ncdf_varid(file_id,'longitude_wgs84')
 ncdf_varget, file_id, ivarid, stn_lon

 ; read the elevation
 ivarid = ncdf_varid(file_id,'elevation')
 ncdf_varget, file_id, ivarid, stn_elev

 ; get time units
 ivar_id = ncdf_varid(file_id,'time')
 ncdf_attget, file_id, ivar_id, 'units', bunits
 cunits = string(bunits)

 ; get the base julian day
 tunit_words = strsplit(string(cunits),' ',/extract)
 tunit_idate = fix(strsplit(tunit_words[2],'-',/extract))
 tunit_ihour = fix(strsplit(tunit_words[3],':',/extract))
 bjulian_stn = julday(tunit_idate[1],tunit_idate[2],tunit_idate[0],tunit_ihour[0],tunit_ihour[1],tunit_ihour[2])

 ; read the time vector
 if(strtrim(tunit_words[0],2) eq 'seconds') then aoff=1.d/86400.d else stop, 'unknown time units'  ; get offset in days
 ncdf_varget, file_id, ivar_id, atime
 djulian_stn = bjulian_stn + atime*aoff

; close file
ncdf_close, file_id

; *****
; (3) DEFINE PRELIMINARY DATA FOR THE INTERPOLATION...
; ****************************************************

; get the number of stations
nSta = n_elements(stn_elev)

; get the size of the grid
nGrid = size(grid_elev, /dimensions)

; get a mask where desire information
ixMask = where(grid_mask eq 1, nMask, complement=jxMask)
if(nMask eq 0)then stop, 'expect need to interpolate somewhere'

; define variables from the station data
cVars_stn = ['ppta','tmp3','dpt3','sol','wnd3sa']
nVars_stn = n_elements(cVars_stn)

; define variables from the daymet grid
cVars_grid = ['pptrate','airtemp','dewtemp','swRadDown','windspd']
nVars_grid = n_elements(cVars_grid)
if(nVars_grid ne nVars_stn)then stop, 'expect same number of variables for the stn and grid'

; define start time
iyyy = 1999
im   = 10
id   = 1
ih   = 1
sjulian=julday(im,id,iyyy,ih,0,0.d)

; identify the start index
xMin = min(abs(djulian_stn - sjulian), iStart_stn)
if(abs(xMin) gt 0.0001d)then stop, 'unable to identify start time'

; get time since the reference time (seconds)
timeSinceRefTimeGrid = (sJulian - bjulian_grid)*86400.d

; define the minimum number of stations to compute the regression
minRegress = 5

; define a very small value (used in the log transform)
verySmall = 0.0000001d

; define the time step (seconds)
dt = 3600.d

; define the number of time steps per day
nPerDay = long(floor(86400.d / dt + 0.5d))

; define the number of days
nDays = long(365*9)

; define the number of time steps
nTime = nDays*nPerDay

; *****
; (4) LOOP THROUGH TIME AND VARIABLES...
; **************************************

; loop through time
for itime=0,nTime-1 do begin

 ; identify time index for the station data
 jTime = itime + iStart_stn

 ; define time variable (seconds
 xTime = double(itime)*dt + timeSinceRefTimeGrid

 ; print progress
 caldat, double(itime)/double(nPerDay) + sJulian, im, id, iyyy, ih
 print, iyyy, im, id, ih, itime, ntime, format='(i4,1x,3(i2,1x),2(i10,1x))'

 ; write time variable
 file_id = ncdf_open(grid_file, /write)
  ivarid = ncdf_varid(file_id,'time')
  ncdf_varput, file_id, ivarid, xTime, offset=[itime], count=[1]
 ncdf_close, file_id

 ; loop through variables
 for ivar=0,nVars_stn-1 do begin

  ; *****
  ; (5) GET DATA SUBSET...
  ; **********************

  ; get the station data
  file_id = ncdf_open(stn_file, /nowrite) 
   ivarid = ncdf_varid(file_id,cVars_stn[iVar])
   ncdf_varget, file_id, ivarid, stnData, offset=[0,jTime], count=[nSta,1]
  ncdf_close, file_id

  ; set the first two stations to missing (don't want to use them)
  stnData[0:1] = -9999.

  ; identify a subset of station data
  iValid = where(stnData gt -100.d, nValid)
  if(nValid lt 2)then stop, 'expect at least 2 valid stations'

  ; get station subset
  dataSubset = stnData[iValid]
  xLatSubset = stn_lat[iValid]
  xLonSubset = stn_lon[iValid]
  elevSubset = stn_elev[iValid]

  ; get the temporary array for the predictions
  xPred = fltarr(nGrid[0],nGrid[1])

  ; get the grid array
  xGrid = fltarr(nGrid[0],nGrid[1])

  ; skip variables where constant
  case cVars_stn[iVar] of

   ; precip
   'ppta': begin
    if(total(dataSubset) lt 0.1d)then begin
     xGrid[ixMask] = mean(dataSubset)
     xGrid[jxMask] = -9999.d  ; outside basin
     iskip=1
    endif else begin
     iskip=0
    endelse
   end

   ; solar radiation
   'sol': begin
    if(total(abs(dataSubset)) lt 0.1d or stddev(dataSubset) lt 0.1d)then begin
     xGrid[ixMask] = mean(dataSubset)
     xGrid[jxMask] = -9999.d ; outside basin
     iskip=1
    endif else begin
     iskip=0
    endelse
   end

   ; everything else
   else: iskip=0

  endcase

  ; *****
  ; (6) INTERPOLATE...
  ; ******************

  ; check if there is a need to interpolate
  if(iskip eq 0)then begin

   ; transform the data (get more normal, and avoid negative values)
   if(cVars_stn[iVar] eq 'ppta' or cVars_stn[iVar] eq 'wnd3sa')then begin
    dataTransform = alog(dataSubset + verySmall)
   endif else begin
    dataTransform = dataSubset
   endelse

   ; compute the fit to the data
   if(nValid ge minRegress)then begin

    ; regress data against elevation
    b = regress(elevSubset, dataTransform, const=a, yfit=stnPred)

    ; interpolate to the grid
    xPred[ixMask] = a + b[0]*grid_elev[ixMask] 

   ; not enough data for a regression -- use the mean
   endif else begin
    stnMean       = mean(dataTransform)
    xPred[ixMask] = stnMean
    stnPred       = replicate(stnMean,nValid)
   endelse

   ; define residuals
   xRes = dataTransform - reform(stnPred)

   ; get the distance array
   xDist = fltarr(nValid)

   ; loop through the grid cells
   for iGrid=0,nGrid[0]-1 do begin
    for jGrid=0,nGrid[1]-1 do begin
 
     ; check that we are within the basin
     if(grid_mask[iGrid,jGrid] eq 1)then begin

      ; get the distance to each station (m)
      for iSta=0,nValid-1 do begin
       xDist[iSta] = map_2points(grid_lon[iGrid,jGrid], grid_lat[iGrid,jGrid], xLonSubset[ista], xLatSubset[iSta], /meters)
      endfor

      ; get the weights
      xTemp = 1.d/xDist
      xWght = xTemp / total(xTemp)

      ; get the grid value (in transformed space)
      xTrans = total(xWght[*] * xRes[*]) + xPred[iGrid,jGrid]      ; detrended inverse distance
      ;xTrans = total(xWght[*] * dataTransform[*])                    ; straight inverse distance

      ; back-transform the data (get more normal, and avoid negative values)
      if(cVars_stn[iVar] eq 'ppta' or cVars_stn[iVar] eq 'wnd3sa')then begin
       xGrid[iGrid,jGrid] = exp(xTrans)
      endif else begin
       xGrid[iGrid,jGrid] = xTrans
      endelse

      ; test
      ;print, 'iGrid, jGrid = ', iGrid, jGrid
      ;print, 'xLon, xLat   = ', grid_lon[iGrid,jGrid], grid_lat[iGrid,jGrid], format='(a,1x,2(f12.3,1x))'
      ;print, 'xDist        = ', xDist,      format='(a,1x,20(f12.3,1x))'
      ;print, 'xWght        = ', xWght,      format='(a,1x,20(f12.3,1x))'

     ; if grid is outside the basin
     endif else begin
      xGrid[iGrid,jGrid] = -9999.d
     endelse   ; grid is outside the basin

    endfor ; looping through the grid elements
   endfor ; looping through the grid elements

   ; test
   ;print, 'xRes          = ', xRes,          format='(a,1x,20(f12.3,1x))'
   ;print, 'dataSubset    = ', dataSubset,    format='(a,1x,20(f12.3,1x))'
   ;print, 'dataTransform = ', dataTransform, format='(a,1x,20(f12.3,1x))'
   ;print, xGrid, format='(20(f9.3,1x))'
   ;stop

  endif  ; if not skipping the grid

  ; *****
  ; (7) WRITE INTERPOLATED GRIDS TO THE NETCDF FILE...
  ; **************************************************

  ; open file for writing
  file_id = ncdf_open(grid_file, /write)

   ; write grid
   ivarid = ncdf_varid(file_id,cVars_grid[iVar])
   ncdf_varput, file_id, ivarid, xGrid, offset=[0,0,itime], count=[nGrid[0],nGrid[0],1]

   ; write the number of valid stations
   ivarid = ncdf_varid(file_id,cVars_grid[iVar]+'_nValid')
   ncdf_varput, file_id, ivarid, nValid, offset=[itime], count=[1]

  ; close file
  ncdf_close, file_id

 endfor  ; (looping through variables)

endfor  ; (looping through time)


stop
end
