pro stn2grid_tollgateMeta

; used to define the grid for the tollgate interpolation

; *****
; (1) DEFINE THE DAYMET GRID...
; *****************************

; define the file path
daymet_path = '/home/mclark/summa/input/tollgate/DayMet_Tile_11912/'

; define the filename (just an example file -- only need the coordinates)
daymet_name = daymet_path + 'dayl_1980.nc'

; open netcdf file for reading
file_id = ncdf_open(daymet_name, /nowrite)

 ; read the x coordinate
 ivarid = ncdf_varid(file_id,'x')
 ncdf_varget, file_id, ivarid, tile_x

 ; read the y coordinate
 ivarid = ncdf_varid(file_id,'y')
 ncdf_varget, file_id, ivarid, tile_y

 ; read the latitude
 ivarid = ncdf_varid(file_id,'lat')
 ncdf_varget, file_id, ivarid, tile_lat

 ; read the longitude
 ivarid = ncdf_varid(file_id,'lon')
 ncdf_varget, file_id, ivarid, tile_lon

; close the file
ncdf_close, file_id


; *****
; (2) GET DAYMET ELEVATION...
; ***************************

; define the header
cHead=''
nHead=6

; define the file
filenm = '/home/mclark/summa/ancillary_data/tollgate/daymet_avgelev_from_nhdplus_lcc.txt'

; open file for reading
openr, in_unit, filenm, /get_lun

 ; loop through header lines
 for iHead=0,nHead-1 do begin
  ; read header
  readf, in_unit, cHead
  cData = strsplit(cHead,' ',/extract)
  ; extract grid info
  case cData[0] of
   'ncols':        nCols     = long(cData[1])
   'nrows':        nRows     = long(cData[1])
   'xllcorner':    xll       = float(cData[1])
   'yllcorner':    yll       = float(cData[1])
   'cellsize':     cSize     = float(cData[1])
   'NODATA_value': ixMissing = long(cData[1])
   else: stop, 'unable to find header value'
  endcase
 endfor  ; end looping through header

 ; extract grid
 daymetElev = fltarr(nCols,nRows)
 readf, in_unit, daymetElev
 ;daymetElev = temporary(reverse(daymetElev,2))

; close file
free_lun, in_unit

; define x and y coordinates
daymet_x = xll + dindgen(nCols)*cSize + cSize/2.d
daymet_y = yll + dindgen(nRows)*cSize + cSize/2.d


; *****
; (3) IDENTIFY DESIRED GRID POINTS...
; ***********************************

; define number of values in the subset
nx_subset = 20
ny_subset = 20

; define the start index
i1 = 110
j1 = 100

; define the end index
i2 = i1 + nx_subset -1
j2 = j1 + ny_subset -1

; define the mask file
imask = intarr(nCols,nRows)
imask[*,*] = 0

; define the file path to the correspondence file
cros_path = '/home/mclark/summa/ancillary_data/tollgate/'

; define the filename (just an example file -- only need the coordinates)
cros_name = cros_path + 'Correspondence.nc'

; open netcdf file for reading
file_id = ncdf_open(cros_name, /nowrite)

 ; read the number of overlaps
 ivarid = ncdf_varid(file_id,'overlaps')
 ncdf_varget, file_id, ivarid, nOverlap

 ; read the i-index
 ivarid = ncdf_varid(file_id,'i_index')
 ncdf_varget, file_id, ivarid, i_index

 ; read the j-index
 ivarid = ncdf_varid(file_id,'j_index')
 ncdf_varget, file_id, ivarid, j_index

; close the NetCDF file
ncdf_close, file_id

; get file dimensions
cros_dims = size(i_index, /dimensions)

; loop through file and identify the i and j index
for jx=0,cros_dims[1] -1 do begin
 for ix=0,nOverlap[jx]-1 do begin
  imask[i_index[ix,jx],j_index[ix,jx]] = 1
 endfor
endfor


; *****
; (4) DEFINE THE NEW NETCDF FILE...
; *********************************

; define the file path
grid_path = '/home/mclark/summa/input/tollgate/'

; define the filename (just an example file -- only need the coordinates)
grid_name = grid_path + 'stn2grid_tollgate.nc'

; define variable names
varnames = [$
           'pptrate',    $  ; precipitation rate (mm/hour)
           'airtemp',    $  ; air temperature (degrees C)
           'dewtemp',    $  ; dewpoint temperature (degrees C)
           'swRadDown',  $  ; incoming shortwave radiation flux (W m-2)'
           'windspd'     ]  ; wind speed (m/s)

; define variable descriptions
var_desc = [$
           'precipitation rate', $
           'air temperature',    $
           'dewpoint temperature', $
           'incoming shortwave radiation flux', $
           'wind speed']

; define the units
var_unit = [$
            'mm h-1', $
            'degrees C', $
            'degrees C', $
            'W m-2',     $
            'm s-1']

; define the base julian day
iy_start = 1983
im_start = 10
id_start = 1
bjulian = julday(im_start,id_start,iy_start,0,0,0.d)
unittxt = 'seconds since '+strtrim(iy_start,2)+'-'+strtrim(im_start,2)+'-'+strtrim(id_start,2)+' 0:0:0.0 -0:00'

; define file
file_id = ncdf_create(strtrim(grid_name,2), /clobber)

 ; define the x and y dimensions
 x_id = ncdf_dimdef(file_id, 'x', nx_subset)
 y_id = ncdf_dimdef(file_id, 'y', ny_subset)

 ; define time dimension
 time_id = ncdf_dimdef(file_id, 'time', /unlimited)

 ; define the x coordinate
 ivarid = ncdf_vardef(file_id, 'x', [x_id], /short)
 ncdf_attput, file_id, ivarid, 'long_name', 'x coordinate'

 ; define the y coordinate
 ivarid = ncdf_vardef(file_id, 'y', [y_id], /short)
 ncdf_attput, file_id, ivarid, 'long_name', 'y coordinate'

 ; define the time variable
 ivarid = ncdf_vardef(file_id, 'time', time_id, /double)
 ncdf_attput, file_id, ivarid, 'units', unittxt, /char

 ; define the mask
 ivarid = ncdf_vardef(file_id, 'mask', [x_id,y_id], /short)
 ncdf_attput, file_id, ivarid, 'long_name', 'tollgate mask'

 ; define the latitude
 ivarid = ncdf_vardef(file_id, 'lat', [x_id,y_id], /float)
 ncdf_attput, file_id, ivarid, 'long_name', 'latitude'
 ncdf_attput, file_id, ivarid, 'units', 'degrees north'

 ; define the longitude
 ivarid = ncdf_vardef(file_id, 'lon', [x_id,y_id], /float)
 ncdf_attput, file_id, ivarid, 'long_name', 'longitude'
 ncdf_attput, file_id, ivarid, 'units', 'degrees east'

 ; define the elevation
 ivarid = ncdf_vardef(file_id, 'elev', [x_id,y_id], /float)
 ncdf_attput, file_id, ivarid, 'long_name', 'elevation'
 ncdf_attput, file_id, ivarid, 'units', 'm'

 ; define data variables
 for ivar=0,n_elements(varnames)-1 do begin
  ivarid = ncdf_vardef(file_id, varnames[ivar], [x_id, y_id, time_id], /float)
  ncdf_attput, file_id, ivarid, 'long_name', strtrim(var_desc[ivar],2), /char
  ncdf_attput, file_id, ivarid, 'units', strtrim(var_unit[ivar],2), /char
  ncdf_attput, file_id, ivarid, '_FillValue', -9999., /float
 endfor

 ; define number of valid stations
 for ivar=0,n_elements(varnames)-1 do begin
  ivarid = ncdf_vardef(file_id, varnames[ivar]+'_nValid', [time_id], /short)
  ncdf_attput, file_id, ivarid, 'long_name', 'number of valid stations', /char
 endfor

 ; end control
 ncdf_control, file_id, /endef

; close the netcdf file
ncdf_close, file_id

; *****
; (5) WRITE METADATA...
; *********************

; open netcdf file for writing
file_id = ncdf_open(grid_name, /write)

 ; write x
 ivarid = ncdf_varid(file_id,'x')
 ncdf_varput, file_id, ivarid, indgen(nx_subset)+i1

 ; write y
 ivarid = ncdf_varid(file_id,'y')
 ncdf_varput, file_id, ivarid, indgen(ny_subset)+j1

 ; write lat
 ivarid = ncdf_varid(file_id,'lat')
 ncdf_varput, file_id, ivarid, tile_lat[i1:i2,j1:j2]

 ; write lon
 ivarid = ncdf_varid(file_id,'lon')
 ncdf_varput, file_id, ivarid, tile_lon[i1:i2,j1:j2]

 ; write elev
 ivarid = ncdf_varid(file_id,'elev')
 ncdf_varput, file_id, ivarid, daymetElev[i1:i2,j1:j2]

 ; write mask
 ivarid = ncdf_varid(file_id,'mask')
 ncdf_varput, file_id, ivarid, imask[i1:i2,j1:j2]

; close the NetCDF file
ncdf_close, file_id




stop
end
