pro forcing_RME

; *****
; (0) DEFINE THE NETCDF FILE...
; *****************************

; define netcdf filename
filename_netcdf = '/home/mclark/summa/input/tollgate/RME_forcing.nc'

; define variable names
varnames = [$
           'pptrate',    $  ; precipitation rate (kg m-2 s-1)
           'SWRadAtm',   $  ; downwelling shortwave radiaiton (W m-2)
           'LWRadAtm',   $  ; downwelling longwave radiation  (W m-2)
           'airtemp',    $  ; air temperature (K)
           'windspd',    $  ; wind speed (m/s)
           'airpres',    $  ; pressure (Pa)
           'spechum']       ; specific humidity (g/g)

; define variable descriptions
var_desc = [$
           'precipitation rate',              $
           'downwelling shortwave radiaiton', $
           'downwelling longwave radiation',  $
           'air temperature',                 $
           'wind speed',                      $
           'air pressure',                    $
           'specific humidity']

; define the units
var_unit = [$
           'kg m-2 s-1',       $
           'W m-2',            $
           'W m-2',            $
           'K',                $
           'm s-1',            $
           'Pa',               $
           'g g-1']

; define station names
stn_names = ['met_sheltered','met_exposed']
nSta = n_elements(stn_names)

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

 ; define other variables
 for ivar=0,n_elements(varnames)-1 do begin
  ivarid = ncdf_vardef(file_id, varnames[ivar], [stn_id, time_id], /float)
  ncdf_attput, file_id, ivarid, 'long_name', strtrim(var_desc[ivar],2), /char
  ncdf_attput, file_id, ivarid, 'units', strtrim(var_unit[ivar],2), /char
  ncdf_attput, file_id, ivarid, '_FillValue', -9999., /float
 endfor

 ; end control
 ncdf_control, file_id, /endef

; close the netcdf file
ncdf_close, file_id

; *****
; (0) READ/WRITE ASCII DATA...
; ****************************

; open netcdf file for writing
file_id = ncdf_open(filename_netcdf, /write)

; define constants
Tfreeze = 273.16

; define desired variables
cdesire = ['CORR', '%Snow', 'Ta', 'RH', 'DPT', 'Tg30', 'ws', 'Si', 'Ii']
ndesire = n_elements(cdesire)

print, cdesire

; define directory
datapath='/home/mclark/summa/validation_data/reynolds_creek/'

; define site names
sitename=['sheltered','exposed']

; define elevations
elevsite=[2049.,2094.]

; define a line of character data
cline_ppt='.'
cline_met='.'

; loop through sites
for isite=0,1 do begin

  ; populate station name
 ivarid = ncdf_varid(file_id,'station')
 ncdf_varput, file_id, ivarid, strtrim(stn_names[isite],2), offset=[0,isite], count=[strlen(strtrim(stn_names[isite],2)),1]

 ; define pressure (constant, Pa)
 apres = 101325. * ( (293.-0.0065*elevsite[isite]) / 293. )^5.256

 ; get precipitation and met filenames
 ppt_filename=datapath + 'ppt_' + sitename[isite] + '.txt'
 met_filename=datapath + 'met_' + sitename[isite] + '.txt'

 ; get number of lines in the file
 nlines_ppt = file_lines(ppt_filename)-1
 nlines_met = file_lines(ppt_filename)-1
 if(nlines_ppt ne nlines_met)then stop, 'number of lines do not match'

 ; define data array
 idate = intarr(5,nlines_ppt)
 adata = fltarr(ndesire,nlines_ppt)

 ; define the data path
 data_path = '/home/mclark/summa/input/tollgate/'

 ; open output file for writing
 openw, out_unit, data_path + 'forcing_' + sitename[isite] + '.txt', /get_lun

 ; open precipitation file for reading
 openr, ppt_unit, ppt_filename, /get_lun

 ; open meteorological file for reading
 openr, met_unit, met_filename, /get_lun

 ; read headers
 readf, ppt_unit, cline_ppt
 readf, met_unit, cline_met

 ; extract headers
 chead_ppt = strsplit(cline_ppt, string(9b), /extract)
 chead_met = strsplit(cline_met, string(9b), /extract)

 ; get number of elements
 nppt = n_elements(chead_ppt)
 nmet = n_elements(chead_met)

 ; define the index in cdesire
 ihead_ppt = intarr(nppt)
 ihead_met = intarr(nmet)

 ; identify the position of the desired data
 for itry=0,nppt-1 do ihead_ppt[itry] = where(strmatch(cdesire,chead_ppt[itry]))
 for itry=0,nmet-1 do ihead_met[itry] = where(strmatch(cdesire,chead_met[itry]))

 ; loop through the lines in the file
 for iline=0,nlines_ppt-1 do begin

  ; read a line of data for precip
  readf, ppt_unit, cline_ppt

  ; extract ppt data into character strings
  cdata_ppt = strsplit(cline_ppt, string(9b), /extract)
  if(n_elements(cdata_ppt) ne nppt)then stop, 'ppt: unexpected number of data elements'

  ; read a line of data for met
  readf, met_unit, cline_met

  ; extract met data into character strings
  cdata_met = strsplit(cline_met, string(9b), /extract)
  if(n_elements(cdata_met) ne nmet)then stop, 'met: unexpected number of data elements'
  
  ; check that the dates match
  for itry=0,4 do if(cdata_met[itry] ne cdata_ppt[itry])then stop, 'dates do not match'

  ; get date
  idate[0:4,iline] = long(cdata_met[0:4])

  ; put data in desired position
  for ippt=0,nppt-1 do if(ihead_ppt[ippt] ge 0)then adata[ihead_ppt[ippt],iline] = float(cdata_ppt[ippt])
  for imet=0,nmet-1 do if(ihead_met[imet] ge 0)then adata[ihead_met[imet],iline] = float(cdata_met[imet])

  ; print progress
  if(idate[2,iline] eq 1) then print, idate[*,iline]
  ;print, idate[*,iline], adata[*,iline], format='(i4,1x,3(i2,1x),i4,1x,9(f12.5,1x))'
  ;if (iline gt 300) then stop

  ; convert temperature and dewpoint from C to K
  atemp = adata[2,iline] + Tfreeze
  dewpt = adata[4,iline] + Tfreeze

  ; convert relative humidity to percent
  relhm = adata[3,iline]*100.

  ; compute specific humidity
  sphum = call_function('DEWPT2SPHM', dewpt, apres)

  ; convert precipitation from mm/hr to kg m-2 s-1
  aprcp = adata[0,iline]/3600.

  ; print processed data to the output file
  printf, out_unit, idate[4,iline], idate[1:3,iline], 0, 0.d, aprcp, adata[7:8,iline], atemp, adata[6,iline], apres, sphum, $
                     adata[1,iline], dewpt, adata[5,iline]+Tfreeze, idate[0,iline], $
                     format='(i4,1x,4(i2,1x),f6.1,1x,e14.4,1x,5(f10.3,1x),e12.3,1x,3(f10.3,1x),i4)'


  ; extract the time information
  iy = idate[4,iline]
  im = idate[1,iline]
  id = idate[2,iline]
  ih = idate[3,iline]

  ; get the julian day
  djulian = julday(im,id,iy,ih,0,0.d)

  ; get the time index
  ix_time = floor((djulian - bjulian)*24.d + 0.5d) - 1L
  ;print, stn_names[ista], ix_time, iy, im, id, ih, format='(a30,1x,i12,1x,i4,1x,3(i2,1x))'

  ; write the time variable
  ivarid = ncdf_varid(file_id, 'time')
  ncdf_varput, file_id, ivarid, double(ix_time+1)*3600.d, offset=ix_time, count=1

  ; get the vector of data
  vData = [aprcp, adata[7:8,iline], atemp, adata[6,iline], apres, sphum]

  ; write data
  for ivar=0,n_elements(varnames)-1 do begin
   ivarid = ncdf_varid(file_id,strtrim(varnames[ivar],2))
   ncdf_varput, file_id, ivarid, vData[ivar], offset=[isite,ix_time], count=[1,1]
  endfor ; looping through variables

 endfor  ; (looping through lines of data)

 free_lun, ppt_unit
 free_lun, met_unit
 free_lun, out_unit

endfor ; looping through sites

; close netcdf file
ncdf_close, file_id

stop
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








