pro basin2ascii_tollgate

; used to make hourly forcing data for all sub-basins above tollgate

; define variables for the basin
cVars_bas = ['pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum']
nVars_bas = n_elements(cVars_bas)

; define the data path
data_path = '/home/mclark/summa/input/tollgate/'

; define basin filename
filenm_bas = data_path + 'grid2basin_tollgate.nc'

; open files
nc_bas = ncdf_open(filenm_bas, /nowrite)

 ; get the basin vector
 ivar_id = ncdf_varid(nc_bas,'basin')
 ncdf_varget, nc_bas, ivar_id, ixHRU
 nBasins = n_elements(ixHRU)

 ; get the time vector
 djulian_bas = call_function('exTime',nc_bas)
 nTime = n_elements(djulian_bas)

 ; define array for all variables
 datArray = dblarr(nVars_bas,nTime)

 ; get the year/month/day of the station data
 ; NOTE small offset to ensure hour 0 is the first day
 caldat, djulian_bas+0.00001d, jm, jd, jyyy, jh, jmin

 ; loop through basins
 for iBas=0,nBasins-1 do begin

  ; define ASCII file for the basin data
  filenm_ASCII = data_path + 'tollgateASCII/basinForcing.' + strtrim(string(ixHRU[iBas],format='(i4.4)'),2) + '.txt'
  print, 'filenm_ASCII = ', filenm_ASCII

  ; open file for writing
  openw, out_unit, filenm_ASCII, /get_lun

  ; loop through variables
  for iVar=0,nVars_bas-1 do begin

   ; get data for given variable
   ivar_id = ncdf_varid(nc_bas,cVars_bas[iVar])
   ncdf_varget, nc_bas, ivar_id, xData, offset=[ibas,0], count=[1,ntime]

   ; populate data vector
   datArray[iVar,*] = xData[0,*]

  endfor  ; looping through variables

  ; write data
  for iTime=0,nTime-1 do begin
   ;print,            jyyy[iTime], jm[iTime], jd[iTime], jh[iTime], jmin[iTime], 0.d, datArray[*,iTime], format='(i4,1x,4(i2,1x),f6.1,1x,e14.4,1x,5(f10.3,1x),e12.3)'
   printf, out_unit, jyyy[iTime], jm[iTime], jd[iTime], jh[iTime], jmin[iTime], 0.d, datArray[*,iTime], format='(i4,1x,4(i2,1x),f6.1,1x,e14.4,1x,5(f10.3,1x),e12.3)'
  endfor

  ; free up the file unit
  free_lun, out_unit

 endfor  ; looping through basins

; close the NetCDF file
ncdf_close, nc_bas

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
