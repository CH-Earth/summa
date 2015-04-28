pro forcing_tollgate

; used to convert ASCII station data to a NetCDF file

; *****
; (1) READ IN STATION METADATA...
; *******************************

; define file
filename = '/home/mclark/summa/input/tollgate/stationData/tollgate_station_meta.csv'

; define line of data
cLine = ''

; define number of stations
nStations = file_lines(filename)-1  ; -1 because of the header

; define variables
keyname = strarr(nStations)
cSite   = strarr(nStations)
stnName = strarr(nStations)
xCoord  = dblarr(nStations)
yCoord  = dblarr(nStations)
zElev   = dblarr(nStations)
aLat    = dblarr(nStations)
aLon    = dblarr(nStations)
cPOR    = strarr(nStations)

; open file for reading
openr, in_unit, filename, /get_lun

 ; read header
 readf, in_unit, cLine
 cHead = strsplit(cLine,',',/extract,count=nVars)

 ; loop through stations
 for iStation=0,nStations-1 do begin

  ; read metadata
  readf, in_unit, cLine
  cData = strsplit(cLine,',',/extract,count=nData)
  if(nData ne nVars)then stop, 'expect nData=nVars'

  ; loop through variables
  for iVar=0,nVars-1 do begin

   ; select variable
   case cHead[iVar] of

    ; process data
    'keyname':          keyname[iStation] = strmid(cData[iVar],0,strpos(cData[iVar],'_'))
    'site':             cSite[iStation]   = cData[iVar]
    'name':             stnName[iStation] = cData[iVar]
    'easting_nad27':    xCoord[iStation]  = double(cData[iVar])
    'northing_nad27':   yCoord[iStation]  = double(cData[iVar])
    'msl':              zElev[iStation]   = double(cData[iVar])
    'latitude_wgs84':   aLat[iStation]    = double(cData[iVar])
    'longitude_wgs84':  aLon[iStation]    = double(cData[iVar])
    'por':              cPOR[iStation]    = cData[iVar]

   endcase
  endfor  ; looping through variables

  ; print results
  print, keyname[iStation], stnName[iStation], aLat[iStation], aLon[iStation], zElev[iStation], format='(a20,1x,a90,1x,3(f9.3,1x))'

 endfor  ; looping through stations

; free up file unit
free_lun, in_unit

; *****
; (2) DEFINE NETCDF FILE...
; *************************

; define the netcdf filename
filename_netcdf = '/home/mclark/summa/input/tollgate/stationData/tollgate_forcing.nc'

; define variable names
varname = ['ppts',$
           'pptu',$
           'ppta',$
           'tmp3',$
           'hum3',$
           'vap3',$
           'dpt3',$
           'sol',$
           'wnd3sa',$
           'wnd3d']

; define long names
vardesc = ['Precipitation, Shielded Raingage', $
           'Precipitation, Unshielded Raingage', $
           'Precipitation, Hamon 1971 Dual Gage Wind Corrected', $
           'Air Temperature', $
           'Relative Humidity', $
           'Vapor Pressure', $
           'Dewpoint Temperature', $
           'Incoming Shortwave Solar Radiation', $
           'Average Wind Speed', $
           'Wind Direction']

; define variable units
varunit = ['mm/hour', $
           'mm/hour', $
           'mm/hour', $
           'degrees C', $
           'decimal percent', $
           'Pa', $
           'degrees C', $
           'W/m2', $
           'm/s', $
           'degrees from N']

; define the base julian day
iy_start = 1983
im_start = 10
id_start = 1
bjulian = julday(im_start,id_start,iy_start,0,0,0.d)
unittxt = 'seconds since '+strtrim(iy_start,2)+'-'+strtrim(im_start,2)+'-'+strtrim(id_start,2)+' 0:0:0.0 -0:00'

; define file
file_id = ncdf_create(strtrim(filename_netcdf,2), /clobber)

 ; define length of string
 str_id = ncdf_dimdef(file_id, 'stringLength', 90)

 ; define station dimension
 stn_id  = ncdf_dimdef(file_id, 'station', nStations)

 ; define time dimension
 time_id = ncdf_dimdef(file_id, 'time', /unlimited)

 ; define the station identifier
 ivarid = ncdf_vardef(file_id, 'station_key', [str_id,stn_id], /char)
 ncdf_attput, file_id, ivarid, 'long_name', 'station key'

 ; define the station name
 ivarid = ncdf_vardef(file_id, 'station_name', [str_id,stn_id], /char)
 ncdf_attput, file_id, ivarid, 'long_name', 'station name'

 ; define the station latitude
 ivarid = ncdf_vardef(file_id, 'latitude_wgs84', [stn_id], /float)
 ncdf_attput, file_id, ivarid, 'long_name', 'station latitude'

 ; define the station latitude
 ivarid = ncdf_vardef(file_id, 'longitude_wgs84', [stn_id], /float)
 ncdf_attput, file_id, ivarid, 'long_name', 'station longitude'

 ; define the station latitude
 ivarid = ncdf_vardef(file_id, 'elevation', [stn_id], /float)
 ncdf_attput, file_id, ivarid, 'long_name', 'station elevation'

 ; define the time variable
 ivarid = ncdf_vardef(file_id, 'time', time_id, /double)
 ncdf_attput, file_id, ivarid, 'units', unittxt, /char

 ; define other variables
 for ivar=0,n_elements(varname)-1 do begin
  ivarid = ncdf_vardef(file_id, varname[ivar], [stn_id, time_id], /float)
  ncdf_attput, file_id, ivarid, 'long_name', strtrim(vardesc[ivar],2), /char
  ncdf_attput, file_id, ivarid, 'units', strtrim(varunit[ivar],2), /char
  ncdf_attput, file_id, ivarid, '_FillValue', -9999., /float
 endfor

 ; end control
 ncdf_control, file_id, /endef

; close the netcdf file
ncdf_close, file_id

; *****
; (3) WRITE METADATA...
; *********************

; open netcdf file for writing
file_id = ncdf_open(filename_netcdf, /write)

; loop through sites
for iStation=0,nStations-1 do begin

 ; populate station key
 ivarid = ncdf_varid(file_id,'station_key')
 ncdf_varput, file_id, ivarid, strtrim(keyname[iStation],2), offset=[0,iStation], count=[strlen(strtrim(keyname[iStation],2)),1]

 ; populate station name
 ivarid = ncdf_varid(file_id,'station_name')
 ncdf_varput, file_id, ivarid, strtrim(stnName[iStation],2), offset=[0,iStation], count=[strlen(strtrim(stnName[iStation],2)),1]

 ; populate latitude
 ivarid = ncdf_varid(file_id,'latitude_wgs84')
 ncdf_varput, file_id, ivarid, aLat[iStation], offset=[iStation], count=[1]

 ; populate longitude
 ivarid = ncdf_varid(file_id,'longitude_wgs84')
 ncdf_varput, file_id, ivarid, aLon[iStation], offset=[iStation], count=[1]

 ; populate elevation
 ivarid = ncdf_varid(file_id,'elevation')
 ncdf_varput, file_id, ivarid, zElev[iStation], offset=[iStation], count=[1]

endfor  ; looping through stations

; close netcdf file
ncdf_close, file_id

; *****
; (4) READ/WRITE ASCII DATA...
; ****************************

; define the number of variables in the NetCDF file
nVars = n_elements(varname)

; define variable
cVar = ['met','ppt']

; define file path
filepath_ascii = '/home/mclark/summa/input/tollgate/stationData/'

; open netcdf file for writing
file_id = ncdf_open(filename_netcdf, /write)

; loop through stations
for iStation=0,nStations-1 do begin

 ; loop through variables
 for iVar=0,1 do begin

  ; define station name
  filename_ascii = filepath_ascii + cVar[iVar] + '/' + strtrim(keyname[iStation],2) + '_' + cVar[iVar] + '.dat'
 
  ; open file for reading
  openr, in_unit, filename_ascii, /get_lun

  ; get number of lines in the file
  nlines = file_lines(filename_ascii)

  ; loop through lines in the file
  for iLine = 0,nLines-1 do begin

   ; read a line of data
   readf, in_unit, cLine

   ; check that the line is a data line
   if(strmid(cLine,0,1) eq '#')then continue

   ; get the header
   if(stregex(cLine, '^[0123456789]',/boolean) eq 0)then begin
    cHead = strsplit(cLine,',',/extract,count=nData)
    continue
   endif

   ; extract the data
   cData = strsplit(cLine,',',/extract,count=nVals)
   if(nVals ne nData)then stop, 'expect nData elements'

   ; initialize precip data
   ppts = -9999.d 
   pptu = -9999.d 
   ppta = -9999.d 

   ; extract met data
   tmp3   = -9999.d 
   hum3   = -9999.d 
   vap3   = -9999.d 
   dpt3   = -9999.d 
   sol    = -9999.d 
   wnd3sa = -9999.d 
   wnd3d  = -9999.d 

   ; process data
   for iData=0,nData-1 do begin

    ; select case
    case cHead[iData] of

     ; do nothing with some elements
     'datetime':
     'wy'      :
     'wd'      :

     ; extract time
     'year'  : iyyy = long(cData[iData])
     'month' : im   = long(cData[iData])
     'day'   : id   = long(cData[iData])
     'hour'  : ih   = long(cData[iData])
     'minute': imin = long(cData[iData])

     ; extract precip data
     'ppts'  : ppts = double(cData[iData])
     'pptu'  : pptu = double(cData[iData])
     'ppta'  : ppta = double(cData[iData])

     ; extract met data
     'tmp3'   : tmp3   = double(cData[iData])
     'hum3'   : hum3   = double(cData[iData])
     'vap3'   : vap3   = double(cData[iData])
     'dpt3'   : dpt3   = double(cData[iData])
     'sol'    : sol    = double(cData[iData])
     'wnd3sa' : wnd3sa = double(cData[iData])
     'wnd3d'  : wnd3d  = double(cData[iData])

     ; check
     else: stop, 'unable to identify data element'

    endcase

   endfor  ; looping through data elements

   ; build the data vector
   xData = [ppts, pptu, ppta, tmp3, hum3, vap3, dpt3, sol, wnd3sa, wnd3d]

   ; get the julian day
   djulian = julday(im,id,iyyy,ih,imin,0.d)
   if(djulian le bjulian)then continue

   ; get the time index
   ix_time = floor((djulian - bjulian)*24.d + 0.5d) - 1L
   print, cVar[iVar], iStation, keyname[iStation], ix_time, iyyy, im, id, ih, imin, format='(a3,1x,i2,1x,a30,1x,i12,1x,i4,1x,4(i2,1x))'

   ; write the time variable
   ivarid = ncdf_varid(file_id, 'time')
   ncdf_varput, file_id, ivarid, double(ix_time+1)*3600.d, offset=ix_time, count=1

   ; write the data vector
   for jVar=0,nVars-1 do begin
    if(xData[jVar] gt -500.d)then begin
     ivarid = ncdf_varid(file_id,varname[jVar])
     ncdf_varput, file_id, ivarid, xData[jVar], offset=[iStation,ix_time], count=[1,1]
    endif
   endfor

  endfor  ; end looping through lines in the file

  ; free up a file unit
  free_lun, in_unit

 endfor  ; end looping through variables
endfor  ; end looping through stations

; close the NetCDF file
ncdf_close, file_id


stop
end
