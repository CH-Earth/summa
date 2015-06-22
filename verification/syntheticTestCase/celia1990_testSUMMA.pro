pro celia1990_testSUMMA

; define plotting parameters
window, 0, xs=1000, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=[0,1,2]

; define the file path
fpath = '/home/mclark/test_summa/summa/output/celia1990/'

; define the file name
fname = 'celia1990_spinup_testSumma.nc'

; define the HRU index
iHRU=0

; define the number of desired time steps
ntime=100

; loop through variables (liquid and matric)
for ivar=0,1 do begin

 ; define variables
 if(ivar eq 0)then cvar='mLayerVolFracLiq'
 if(ivar eq 1)then cvar='mLayerMatricHead'

 ; make a base plot
 if(ivar eq 0) then begin
  plot, indgen(5), xrange=[0,0.6], yrange=[0.1,0.22], xstyle=1, ystyle=1, $
   xtitle='Depth (m)', ytitle='Volumetric liquid water content (-)', /nodata
 endif else begin
  plot, indgen(5), xrange=[0,0.6], yrange=[-10,0], xstyle=1, ystyle=1, $
   xtitle='Depth (m)', ytitle='Pressure head (m)', /nodata
 endelse

 ; open file
 nc_file = ncdf_open(fpath+fname, /nowrite)

 ; extract the time vector
 ivar_id = ncdf_varid(nc_file,'time')
 ncdf_varget, nc_file, ivar_id, atime

 ; get the number of layers
 ivar_id = ncdf_varid(nc_file,'nLayers')
 ncdf_varget, nc_file, ivar_id, nLayers, offset=[iHRU,0], count=[1,ntime]

 ; get the number of snow layers
 ivar_id = ncdf_varid(nc_file,'nSnow')
 ncdf_varget, nc_file, ivar_id, nSnow, offset=[iHRU,0], count=[1,ntime]

 ; get the number of soil layers
 ivar_id = ncdf_varid(nc_file,'nSoil')
 ncdf_varget, nc_file, ivar_id, nSoil, offset=[iHRU,0], count=[1,ntime]

 ; get the start index for mLayers
 ivar_id = ncdf_varid(nc_file,'midSoilStartIndex')
 ncdf_varget, nc_file, ivar_id, midSoilStartIndex, offset=[iHRU,0], count=[1,ntime]

 ; loop through time
 for itime=0,ntime-1 do begin

  ; define color
  icol = fix((float(itime+1)/52.)*200.) + 50

  ; get the mid-point of each layer
  ivar_id = ncdf_varid(nc_file,'mLayerHeight')
  ncdf_varget, nc_file, ivar_id, mLayerHeight, offset=[iHRU,(midSoilStartIndex[0,itime]+nSnow[0,itime])-1], count=[1,nSoil[0,itime]]

  ; get the desired variable for all layers
  ivar_id = ncdf_varid(nc_file,cvar)
  ncdf_varget, nc_file, ivar_id, avar, offset=[iHRU,midSoilStartIndex[0,itime]-1], count=[1,nSoil[0,itime]]

  ; plot the data
  ;if(itime eq 9 or itime eq 27 or itime eq 51)then begin
  if(itime eq 10 or itime eq 32 or itime eq 49)then begin
   oplot, mLayerHeight[0,*], avar[0,*], color=icol
   oplot, mLayerHeight[0,*], avar[0,*], color=icol, psym=sym(1)
   print, 'time = ', long(itime+1)*1800L
   print, reform(avar[0,*])
  endif

 endfor ; loop through time

 ; close NetCDF file
 ncdf_close, nc_file

endfor  ; loop through variables

; make a figure
write_png, 'zFigures/celia1990_testSUMMA.png', tvrd(true=1)

stop
end
