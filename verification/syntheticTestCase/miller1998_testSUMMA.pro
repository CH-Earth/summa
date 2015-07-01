pro miller1998_testSUMMA

; define plotting parameters
window, 0, xs=1000, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=[0,1,2]

; define the HRU
iHRU=0

; define the maximum soil depth
maxDepth = 10.d

; define the time step
dTime = 900.d

; define the desired time (units of seconds)
xTime = [0.18d, 2.25d, 1.d] * 86400.d

; define the soil depth
sDepth = [10.d,5.d,2.d]

; define the path to the graphics file
gpath = '/home/mclark/test_summa/summa/verification/zFigures/'

; define the name of the graphics file
gname = 'syntheticTestCase_miller1998.png'

; define the file path
fpath = '/home/mclark/test_summa/summa/output/syntheticTestCases/miller1998/'

; define experiment names
expName = ['millerSand','millerLoam','millerClay']

; define file suffix
fSuff = '_spinup_testSumma.nc'

; define colors
icol = [80,210,150]

; loop through variables (liquid and matric)
for ivar=0,1 do begin

 ; define variables
 if(ivar eq 0)then cvar='mLayerVolFracLiq'
 if(ivar eq 1)then cvar='mLayerMatricHead'

 ; make a base plot
 if(ivar eq 0) then begin
  plot, indgen(5), yrange=[maxDepth,0], xrange=[0.05,0.45], xstyle=1, ystyle=1, $
   ytitle='Depth (m)', xtitle='Volumetric liquid water content (-)', /nodata
 endif else begin
  plot, indgen(5), yrange=[maxDepth,0], xrange=[-5,0.5], xstyle=1, ystyle=1, $
   ytitle='Depth (m)', xtitle='Pressure head (m)', /nodata
 endelse

 ; loop through files
 for ifile=0,n_elements(expName)-1 do begin
  
  ; define file names
  filenm = expName[ifile] + fSuff

  ; define the desired time step
  nTime = floor(xTime[ifile]/dTime + 0.5d)

  ; open file
  nc_file = ncdf_open(fpath+filenm, /nowrite)

  ; extract the time vector
  ivar_id = ncdf_varid(nc_file,'time')
  ncdf_varget, nc_file, ivar_id, atime

  ; get the number of soil layers
  ivar_id = ncdf_varid(nc_file,'nSoil')
  ncdf_varget, nc_file, ivar_id, nSoil, offset=[iHRU,ntime], count=[1,1]

  ; get the start index for mLayers
  ivar_id = ncdf_varid(nc_file,'midSoilStartIndex')
  ncdf_varget, nc_file, ivar_id, midSoilStartIndex, offset=[iHRU,ntime], count=[1,1]

  ; get the mid-point of each layer
  ivar_id = ncdf_varid(nc_file,'mLayerHeight')
  ncdf_varget, nc_file, ivar_id, mLayerHeight, offset=[iHRU,midSoilStartIndex-1], count=[1,nSoil]

  ; get the desired variable for all layers
  ivar_id = ncdf_varid(nc_file,cvar)
  ncdf_varget, nc_file, ivar_id, avar, offset=[iHRU,midSoilStartIndex-1], count=[1,nSoil]

  ; plot the data
  oplot, avar[0,*], mLayerHeight[0,*] + (maxDepth - sDepth[ifile]), color=icol[ifile], thick=3
  print, 'time = ', long(ntime+1)*long(dTime)
  print, reform(avar[0,*])

  ; close NetCDF file
  ncdf_close, nc_file

 endfor  ; loop through files
endfor  ; loop through variables

; make a figure
write_png, gpath+gname, tvrd(true=1)

stop
end
