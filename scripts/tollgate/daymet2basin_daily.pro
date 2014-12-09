pro daymet2basin_daily

; used to interpolate daymet data (daily time step) to sub-basins

; *****
; (0) GET INFORMATION FROM THE NETWORK TOPOLOGY FILE...
; *****************************************************

; declare network topology file
filenm = '/home/mclark/summa/ancillary_data/tollgate/Reynolds_Network_Topology.nc'
print, filenm

; open file
nc_file = ncdf_open(filenm, /nowrite)

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
; (1) GET THE INFORMATION FROM THE CORRESPONDENCE FILE...
; *******************************************************

; declare network topology file
filenm = '/home/mclark/summa/ancillary_data/tollgate/Correspondence.nc'
print, filenm

; open file
nc_file = ncdf_open(filenm, /nowrite)

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
; (2) DEFINE FILE FOR THE RE-MAPPED DATA...
; *****************************************

; define variables
cDesire = ['prcp','tmax','tmin','srad','dayl','vp']

; define variable names
cVarName = ['Precipitation','Maximum temperature','Minimum temperature','Daylight average shortwave radiation','Day length','Vapor pressure']

; define units
cVarUnits = ['mm/day','degrees C','degrees C','W/m2','s','Pa']

; get the number of basins
nBasins = n_elements(polyid)

; define the NetCDF file
filenm = '/home/mclark/summa/input/tollgate/tollgateBasin_daily.nc'

; create the NetCDF file
nc_file = ncdf_create(filenm, /clobber)

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
 ncdf_attput, nc_file, ivar_id, 'units', 'days since 1980-1-1 0:0:0.0 -0:00'

 ; define data variable for the basin elevation
 ivar_id = ncdf_vardef(nc_file,'basinElev',[ibas_id], /double)
 ncdf_attput, nc_file, ivar_id, 'long_name', 'basin elevation'
 ncdf_attput, nc_file, ivar_id, 'units', 'm'

 ; define data variable for the basin area
 ivar_id = ncdf_vardef(nc_file,'basinArea',[ibas_id], /double)
 ncdf_attput, nc_file, ivar_id, 'long_name', 'basin area'
 ncdf_attput, nc_file, ivar_id, 'units', 'm2'

 ; define data variables
 for ivar=0,n_elements(cDesire)-1 do begin

  ivar_id = ncdf_vardef(nc_file,cDesire[ivar], [ibas_id, itim_id], /double)
  ncdf_attput, nc_file, ivar_id, 'long_name', cVarName[ivar]
  ncdf_attput, nc_file, ivar_id, 'units', cVarUnits[ivar]

 endfor

; close the NetCDF file
ncdf_close, nc_file

; *****
; (3) POPULATE FILE WITH DAYMET DATA...
; *************************************

; define desired date range
iy1 = 1980
iy2 = 2013

; define the path to the daymet files
daymet_path = '/home/mclark/summa/input/tollgate/DayMet_Tile_11912/'

; define number of days in daymet
nDays = 365  ; note: daymet does not have leap years (correct later)

; define the number of leap years
nLeap = 0

; open the NetCDF file for the basin data
nc_basin = ncdf_open(filenm, /write)

 ; loop through the years
 for iyear=iy1,iy2 do begin

  ; identify the leap year
  if( (399+(iyear mod 400))/400 - (3+(iyear mod 4))/4 eq 1 or iyear eq 2000)then iLeap=1 else iLeap=0
  
  ; define time offset for data write
  ixTime = (iyear - iy1)*nDays + nLeap ; increment leap year later
  print, ixTime, iLeap, nLeap

  ; increment leap year
  ; NOTE: this must be done before defining the time offset
  if(iLeap eq 1)then nLeap = nLeap + 1

  ; loop through desired variables
  for iDesire=0,n_elements(cDesire)-1 do begin

   ; define the daymet file
   daymet_file = daymet_path + cDesire[iDesire] + '_' + strtrim(iyear,2) + '.nc'
   print, 'daymet_file = ', daymet_file

   ; open up the NetCDF file for the daymet data
   nc_daymet = ncdf_open(daymet_file, /nowrite)

    ; loop through basins
    for iBas=0,nBasins-1 do begin

     ; define the basin data
     bData = replicate(0.d, nDays) 

     ; check that the weights sum to 1
     sWeight = total(weight[0:nOverlap[ibas]-1,ibas])
     if(abs(1.d - sWeight) gt 0.001d)then stop, 'weights do not sum to one'

     ; loop through the desired grid points
     for ixOverlap=0,nOverlap[iBas]-1 do begin

      ; get the x and y indices for the daymet grid
      ixDaymet = daymet_i[ixOverlap,ibas]-1
      jyDaymet = daymet_j[ixOverlap,ibas]-1
      ;print, ixDaymet, jyDaymet

      ; get the data
      ivar_id = ncdf_varid(nc_daymet,cDesire[iDesire])
      ncdf_varget, nc_daymet, ivar_id, xData, offset=[ixDaymet,jyDaymet,0], count=[1,1,nDays]

      ; add data to bdata
      bData[*] = bData[*] + xData[0,0,*]*weight[ixOverlap,ibas]

     endfor  ; looping through overlapping gridpoints

     ; NOTE: dec 31 discarded in leap years -- use persistence
     if(iLeap eq 1)then bData = temporary([bData[0:nDays-1],bData[nDays-1]])

     ; write basin coordinate
     ivar_id = ncdf_varid(nc_basin,'basin')
     ncdf_varput, nc_basin, ivar_id, polyid[ibas], offset=[ibas], count=[1]

     ; write basin elevation
     ivar_id = ncdf_varid(nc_basin,'basinElev')
     ncdf_varput, nc_basin, ivar_id, basinElev[ibas], offset=[ibas], count=[1]

     ; write basin area
     ivar_id = ncdf_varid(nc_basin,'basinArea')
     ncdf_varput, nc_basin, ivar_id, basinArea[ibas], offset=[ibas], count=[1]

     ; write time coordinate
     ivar_id = ncdf_varid(nc_basin,'time')
     ncdf_varput, nc_basin, ivar_id, indgen(nDays+iLeap)+ixTime, offset=[ixtime], count=[nDays+iLeap]

     ; write data to the output file
     ivar_id = ncdf_varid(nc_basin,cDesire[iDesire])
     ncdf_varput, nc_basin, ivar_id, bData, offset=[ibas,ixTime], count=[1,nDays+iLeap]

    endfor  ; looping through basins

   ; close the Daymet file
   ncdf_close, nc_daymet

  endfor  ; looping through desired variables
 endfor  ; looping through the years

; close the basin NetCDF file
ncdf_close, nc_basin

stop
end
