pro read_groundHeatFlux

; define the HRU
iHRU=0

; define file name
file_name = '/home/mclark/summa/output/plumber/orig.nc'

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

 ; get the ground heat flux
 for iTime=0,nTime_mod-1 do begin

  ; get the vector of measurment height
  ivar_id = ncdf_varid(nc_file,'iLayerHeight')
  ncdf_varget, nc_file, ivar_id, xHeight, offset=[iHRU,ifcTotoStartIndex[iTime]-1], count=[1,nLayers[iTime]+1]

  ; get the vector of the heat flux
  ivar_id = ncdf_varid(nc_file,'iLayerNrgFlux')
  ncdf_varget, nc_file, ivar_id, xFlux, offset=[iHRU,ifcTotoStartIndex[iTime]-1], count=[1,nLayers[iTime]+1]

  ; print progress
  ;print, itime, ifcTotoStartIndex[iTime], xHeight[iHRU,*], format='(2(i6,1x),20(f13.5,1x))'
  print, itime, ifcTotoStartIndex[iTime], xFlux[iHRU,*], format='(2(i6,1x),10(f13.5,1x))'

 endfor  ; looping through time

; close the netcdf file
ncdf_close, nc_file



end
