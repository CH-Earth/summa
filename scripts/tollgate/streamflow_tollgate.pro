pro streamflow_tollgate

; *****
; (0) DEFINE THE NETCDF FILE...
; *****************************

; define netcdf filename
filename_netcdf = '/home/mclark/summa/validation_data/reynolds_creek/tollgate_streamflow.nc'

; define number of stations
nSta = 1

; define the base julian day
iy_start = 1983
im_start = 10
id_start = 1
bjulian = julday(im_start,id_start,iy_start,0,0,0.d)
unittxt = 'seconds since '+strtrim(iy_start,2)+'-'+strtrim(im_start,2)+'-'+strtrim(id_start,2)+' 0:0:0.0 -0:00'

; define file
file_id = ncdf_create(strtrim(filename_netcdf,2), /clobber)

 ; define length of string
 str_id = ncdf_dimdef(file_id, 'stringLength', 30)

 ; define station dimension
 stn_id  = ncdf_dimdef(file_id, 'station', nSta)

 ; define time dimension
 time_id = ncdf_dimdef(file_id, 'time', /unlimited)

 ; define the station variable
 ivarid = ncdf_vardef(file_id, 'station', [str_id,stn_id], /char)
 ncdf_attput, file_id, ivarid, 'long_name', 'station name'

 ; define the time variable
 ivarid = ncdf_vardef(file_id, 'time', time_id, /double)
 ncdf_attput, file_id, ivarid, 'units', unittxt, /char

 ; define streamflow
 ivarid = ncdf_vardef(file_id, 'streamflow', [stn_id, time_id], /float)
 ncdf_attput, file_id, ivarid, 'long_name', 'streamflow', /char
 ncdf_attput, file_id, ivarid, 'units', 'm3/s', /char
 ncdf_attput, file_id, ivarid, '_FillValue', -9999., /float

 ; end control
 ncdf_control, file_id, /endef

; close the netcdf file
ncdf_close, file_id

; *****
; (0) READ/WRITE ASCII DATA...
; ****************************

; open netcdf file for writing
file_id = ncdf_open(filename_netcdf, /write)

; define directory
datapath='/home/mclark/summa/validation_data/reynolds_creek/'

; define site names
sitename=['tollgate']

; define a line of character data
cline='.'

; loop through sites
for isite=0,n_elements(sitename)-1 do begin

 ; populate station name
 ivarid = ncdf_varid(file_id,'station')
 ncdf_varput, file_id, ivarid, strtrim(sitename[isite],2), offset=[0,isite], count=[strlen(strtrim(sitename[isite],2)),1]

 ; open file for reading
 filename = datapath + 'tollgate_flow.dat'
 openr, in_unit, filename, /get_lun

 ; get number of lines in the file
 nlines = file_lines(filename)

 ; loop through lines in the file
 for iLine = 0,nLines-1 do begin

  ; read a line of data
  readf, in_unit, cLine

  ; check that the line is a data line
  if(strmid(cLine,0,1) eq '#')then continue

  ; get the header
  if(stregex(cLine, '[0123456789]',/boolean) eq 0)then begin
   cHead = strsplit(cLine,',',/extract,count=nData)
   continue
  endif

  ; extract the data
  cData = strsplit(cLine,',',/extract,count=nVals)
  if(nVals ne nData)then stop, 'expect nData elements'

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

    ; extract data
    'qcms'  : qcms = double(cData[iData])

    ; check
    else: stop, 'unable to identify data element'

   endcase

  endfor  ; looping through data elements

  ; get the julian day
  djulian = julday(im,id,iyyy,ih,imin,0.d)
  if(djulian le bjulian)then continue

  ; get the time index
  ix_time = floor((djulian - bjulian)*24.d + 0.5d) - 1L
  print, sitename[isite], ix_time, iyyy, im, id, ih, imin, qcms, format='(a30,1x,i12,1x,i4,1x,4(i2,1x),f20.5)'

  ; write the time variable
  ivarid = ncdf_varid(file_id, 'time')
  ncdf_varput, file_id, ivarid, double(ix_time+1)*3600.d, offset=ix_time, count=1

  ; write data
  ivarid = ncdf_varid(file_id,'streamflow')
  ncdf_varput, file_id, ivarid, qcms, offset=[isite,ix_time], count=[1,1]

 endfor  ; (looping through lines of data)

 ; close file unit
 free_lun, in_unit

endfor ; looping through sites

; close netcdf file
ncdf_close, file_id

stop
end
