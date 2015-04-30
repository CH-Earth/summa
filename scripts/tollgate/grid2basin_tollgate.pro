pro grid2basin_tollgate

; used to make hourly forcing data for all sub-basins above tollgate

; *****
; (1) PRELIMINARIES: DEFINE VARIABLES ETC...
; ******************************************

; define the data path
data_path = '/home/mclark/summa/input/tollgate/'

; define gridded filename
filenm_grid = data_path + 'stn2grid_tollgate.nc'

; define station filename
filenm_stn = data_path + 'RME_forcing.nc'

; define the hourly time step
dt_hour = 3600.d

; open files
nc_grid = ncdf_open(filenm_grid, /nowrite)
nc_stn = ncdf_open(filenm_stn, /nowrite)

 ; get the time vectors
 djulian_grid = call_function('exTime',nc_grid)
 djulian_stn = call_function('exTime',nc_stn)

 ; get the daymet x coordinate
 ivar_id = ncdf_varid(nc_grid,'x')
 ncdf_varget, nc_grid, ivar_id, ix_grid

 ; get the daymet y coordinate
 ivar_id = ncdf_varid(nc_grid,'y')
 ncdf_varget, nc_grid, ivar_id, jy_grid

; close netcdf files
ncdf_close, nc_grid
ncdf_close, nc_stn

; *****
; (2) GET INFORMATION FROM THE NETWORK TOPOLOGY FILE...
; *****************************************************

; declare network topology file
filenm_network = '/home/mclark/summa/ancillary_data/tollgate/Reynolds_Network_Topology.nc'
print, filenm_network

; open file
nc_file = ncdf_open(filenm_network, /nowrite)

 ; get the hruid
 ivar_id = ncdf_varid(nc_file,'hruid')
 ncdf_varget, nc_file, ivar_id, ixHRU

 ; get the upstream area
 ivar_id = ncdf_varid(nc_file,'Elev_Avg')
 ncdf_varget, nc_file, ivar_id, Elev_Avg

 ; get the index of the downstream segment
 ivar_id = ncdf_varid(nc_file,'Basin_Area')
 ncdf_varget, nc_file, ivar_id, Basin_Area

; close file
ncdf_close, nc_file

; *****
; (3) GET THE INFORMATION FROM THE CORRESPONDENCE FILE...
; *******************************************************

; declare network topology file
filenm_grid2bas = '/home/mclark/summa/ancillary_data/tollgate/Correspondence.nc'
print, filenm_grid2bas

; open file
nc_file = ncdf_open(filenm_grid2bas, /nowrite)

 ; get the id of each polygon
 ivar_id = ncdf_varid(nc_file,'polyid')
 ncdf_varget, nc_file, ivar_id, polyid

 ; get the number of grid overlaps
 ivar_id = ncdf_varid(nc_file,'overlaps')
 ncdf_varget, nc_file, ivar_id, noverlap

 ; get the i-index
 ivar_id = ncdf_varid(nc_file,'i_index')
 ncdf_varget, nc_file, ivar_id, daymet_i

 ; get the j-index
 ivar_id = ncdf_varid(nc_file,'j_index')
 ncdf_varget, nc_file, ivar_id, daymet_j

 ; get the weight
 ivar_id = ncdf_varid(nc_file,'weight')
 ncdf_varget, nc_file, ivar_id, weight

; close file
ncdf_close, nc_file

; define the number of basins
nBasins = n_elements(polyid)

; define basin elev and area
basinElev = dblarr(nBasins)
basinArea = dblarr(nBasins)

; put area and elevation data into the polygon array
for iRch=0,n_elements(Basin_Area)-1 do begin
 iMatch = where(ixHRU[iRch] eq polyid, nMatch)
 if(nMatch eq 1)then begin
  basinElev[iMatch] = Elev_Avg[iRch]
  basinArea[iMatch] = Basin_Area[iRch]
 endif
endfor

; *****
; (4) DEFINE FILE FOR THE RE-MAPPED DATA...
; *****************************************

; define the NetCDF file
filenm_bas = '/home/mclark/summa/input/tollgate/grid2basin_tollgate.nc'

; define named variables
ixBas_pptrate   = 0
ixBas_SWRadAtm  = 1
ixBas_LWRadAtm  = 2
ixBas_airtemp   = 3
ixBas_windspd   = 4
ixBas_airpres   = 5
ixBas_spechum   = 6

; define variables for the basin
cVars_bas = ['pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum']
nVars_bas = n_elements(cVars_bas)

; define variable names
cVarName = ['Precipitation','Downward shortwave radiation','Downward longwave radiation','Air temperature','Wind speed','Air pressure','Specific humidity']

; define units
cVarUnits = ['kg m-2 s-1','W m-2','W m-2','K','m s-1','Pa','g g-1']

; get the time units from the grid file (use the same units)
nc_grid = ncdf_open(filenm_grid, /nowrite)
 ivar_id = ncdf_varid(nc_grid,'time')
 ncdf_attget, nc_grid, ivar_id, 'units', bunits
 cunits = string(bunits)
ncdf_close, nc_grid

; get the number of basins
nBasins = n_elements(polyid)

; create the NetCDF file
nc_file = ncdf_create(filenm_bas, /clobber)

 ; define the number of basins
 ibas_id = ncdf_dimdef(nc_file,'basin',nBasins)

 ; define the number of time steps
 itim_id = ncdf_dimdef(nc_file,'time',/unlimited)

 ; define coordinate variable for the basin
 ivar_id = ncdf_vardef(nc_file,'basin', [ibas_id], /long)
 ncdf_attput, nc_file, ivar_id, 'long_name', 'basin index'

 ; define coordinate variable for time
 ivar_id = ncdf_vardef(nc_file,'time', [itim_id], /long)
 ncdf_attput, nc_file, ivar_id, 'long_name', 'time since time reference'
 ncdf_attput, nc_file, ivar_id, 'units', strtrim(cunits,2)

 ; define data variable for the basin elevation
 ivar_id = ncdf_vardef(nc_file,'basinElev',[ibas_id], /double)
 ncdf_attput, nc_file, ivar_id, 'long_name', 'basin elevation'
 ncdf_attput, nc_file, ivar_id, 'units', 'm'

 ; define data variable for the basin area
 ivar_id = ncdf_vardef(nc_file,'basinArea',[ibas_id], /double)
 ncdf_attput, nc_file, ivar_id, 'long_name', 'basin area'
 ncdf_attput, nc_file, ivar_id, 'units', 'm2'

 ; define data variables
 for ivar=0,n_elements(cVars_bas)-1 do begin

  ivar_id = ncdf_vardef(nc_file,cVars_bas[ivar], [ibas_id, itim_id], /double)
  ncdf_attput, nc_file, ivar_id, 'long_name', cVarName[ivar]
  ncdf_attput, nc_file, ivar_id, 'units', cVarUnits[ivar]

 endfor

 ; end the file definitions
 ncdf_control, nc_file, /endef

 ; write the basin id
 ivar_id = ncdf_varid(nc_file,'basin')
 ncdf_varput, nc_file, ivar_id, polyid

 ; write the basin elevation
 ivar_id = ncdf_varid(nc_file,'basinElev')
 ncdf_varput, nc_file, ivar_id, basinElev

 ; write the basin area
 ivar_id = ncdf_varid(nc_file,'basinArea')
 ncdf_varput, nc_file, ivar_id, basinArea

; close the NetCDF file
ncdf_close, nc_file

; *****
; (4) LOOP THROUGH TIME...
; ************************

; define named variables
ixGrid_pptrate   = 0
ixGrid_airtemp   = 1
ixGrid_dewtemp   = 2
ixGrid_swRadDown = 3
ixGrid_windspd   = 4

; define desired variables from the grid
cVars_grid = ['pptrate','airtemp','dewtemp','swRadDown','windspd']
nVars_grid = n_elements(cVars_grid)

; get the year/month/day of the station data
; NOTE small offset to ensure hour 0 is the first day
caldat, djulian_stn-0.00001d, jm, jd, jyyy, jh, jmin

; define the index of the desired site
; NOTE: 0 = sheltered site in RME
iStn = 0

; open grid and station files for reading
nc_grid = ncdf_open(filenm_grid, /nowrite)
nc_stn = ncdf_open(filenm_stn, /nowrite)

; open the basin file for writing
nc_bas = ncdf_open(filenm_bas, /write)

; get the time vector from the grid file (use the same time vector)
ivar_id = ncdf_varid(nc_grid,'time')
ncdf_varget, nc_grid, ivar_id, atime

; loop through time
for itime=0,n_elements(djulian_grid)-1 do begin

 ; write time
 ivar_id = ncdf_varid(nc_bas,'time')
 ncdf_varput, nc_bas, ivar_id, atime[itime], offset=[itime], count=[1]

 ; identify the match with the station data
 caldat, djulian_grid[iTime], im, id, iyyy, ih
 iMatch = where(iyyy eq jyyy and im eq jm and id eq jd and ih eq jh, nMatch)
 if(nMatch ne 1)then stop, 'expect a single match'
 print, iyyy, im, id, ih, format='(i4,1x,3(i2,1x))'

 ; get the longwave radiation data from Reynolds Mountain East
 ivar_id = ncdf_varid(nc_stn,'LWRadAtm')
 ncdf_varget, nc_stn, ivar_id, xLWRadAtm, offset=[iStn,iMatch], count=[1,1]

 ; loop through basins
 for iBas=0,nBasins-1 do begin

  ; *****
  ; (5) INTERPOLATE DATA TO BASINS...
  ; *********************************

  ; check that the weights sum to 1
  sWeight = total(weight[0:nOverlap[ibas]-1,ibas])
  if(abs(1.d - sWeight) gt 0.001d)then stop, 'weights do not sum to one'

  ; define the basin data
  bData_orig = replicate(0.d, nVars_grid)
  bData_conv = dblarr(nVars_bas)

  ; loop through desired variables from the grid
  for iVar=0,nVars_grid-1 do begin

   ; loop through the desired grid points
   for ixOverlap=0,nOverlap[iBas]-1 do begin

    ; get the x and y indices for the daymet grid
    ixDaymet = where(daymet_i[ixOverlap,ibas] eq ix_grid, nxMatch)
    jyDaymet = where(daymet_j[ixOverlap,ibas] eq jy_grid, nyMatch)
    if(nxMatch ne 1 or nyMatch ne 1)then stop, 'expect a single match for the i and j indices'

    ; get the data
    ivar_id = ncdf_varid(nc_grid,cVars_grid[iVar])
    ncdf_varget, nc_grid, ivar_id, xData, offset=[ixDaymet[0],jyDaymet[0],itime], count=[1,1,1]
    ;print, cVars_grid[iVar], ixDaymet[0], jyDaymet[0], xData, format='(a15,1x,2(i2,1x),f9.3,1x)'

    ; add data to bdata
    bData_orig[iVar] = bData_orig[iVar] + reform(xData)*weight[ixOverlap,ibas]

   endfor  ; looping through overlapping gridpoints

  endfor  ; looping through variables

  ; *****
  ; (6) GET DESIRED VARIABLES FOR THE MODEL RUNS...
  ; ***********************************************

  ; define precipitation (mm/h --> kg m-2 s-1)
  bData_conv[ixBas_pptrate] = bData_orig[ixGrid_pptrate]/dt_hour 

  ; define shortwave radiation
  bData_conv[ixBas_SWRadAtm] = bData_orig[ixGrid_swRadDown]

  ; copy pver longwave radiation
  bData_conv[ixBas_LWRadAtm] = reform(xLWRadAtm)

  ; define air temperature (degrees C --> K)
  bData_conv[ixBas_airtemp] = bData_orig[ixGrid_airtemp] + 273.16d

  ; define wind speed
  bData_conv[ixBas_windspd] = bData_orig[ixGrid_windspd]

  ; define pressure
  bData_conv[ixBas_airpres] = 101325.d * ( (293.d - 0.0065d * basinElev[ibas]) / 293.d )^5.256d

  ; define specific humidity
  bData_conv[ixBas_spechum] = call_function('DEWPT2SPHM', bData_orig[ixGrid_dewtemp]+273.16d, bData_conv[ixBas_airpres])

  ; print results
  ;print, 'bData_orig = ', bData_orig, format='(a,10(f9.3,1x))'
  ;print, 'bData_conv = ', bData_conv, format='(a,10(f9.3,1x))'

  ; *****
  ; (7) WRITE DATA TO THE NETCDF FILE...
  ; ************************************
 
  ; write variables
  for iVar=0,nVars_bas-1 do begin
   ivar_id = ncdf_varid(nc_bas,cVars_bas[iVar])
   ncdf_varput, nc_bas, ivar_id, bData_conv[iVar], offset=[ibas,itime], count=[1,1]
  endfor

 endfor  ; looping through basins

endfor  ; looping through time

; close grid and station files
ncdf_close, nc_grid
ncdf_close, nc_stn

; close basin file
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
