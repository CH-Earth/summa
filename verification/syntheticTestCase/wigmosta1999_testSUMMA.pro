pro wigmosta1999_testSUMMA

; define plotting parameters
window, 0, xs=1200, ys=700, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,3,0,1]

; define soil parameters
tan_slope    = 0.3d  ; -
TOPMODEL_exp = 1.d   ; -
surfHydCond  = 0.3d  ; m h-1
soilDepth    = 1.5d  ; m
porosity     = 0.25d ; -

; define forcing
qRain = 0.002 ; m/h

; define number of spatial elements for the kinematic predictions
nCoord = 500

; define the number of hillslope elements for the model simulations
nHill = 50

; define the total length of the hillslope
totalLength = 50.d  ; m

; define hillslope width
hillWidth = 100.d ; m

; define the domain area (m2)
domainArea = hillWidth*totalLength

; define the x-coordinate for the kinematic predictions
xCoord = totalLength*(dindgen(nCoord)+0.5d)/double(nCoord)

; define the x-xoordinate for the model simulations
xSlope = totalLength*(dindgen(nHill)+0.5d)/double(nHill)

; define the file path
fpath = '/home/mclark/test_summa/summa/output/wigmosta1999/'

; define the number of HRUs
nHRU=50

; *************************************************************************************************************
; * PART 1: PLOT HILLSLOPE STORAGE AT DIFFERENT OUTPUT TIMES
; *************************************************************************************************************

; define the total storage
totalStorage = dblarr(nHRU)

; define desired output times
;outTime = [35,70,140,400]
outTime = [35,70,140]

; define x and y titles
ytit = [' ', 'Drainable storage (m)', ' ']
xtit = [' ', ' ', 'Downslope distance (m)']
;ytit = [' ', 'Drainable storage (m)', ' ', ' ']
;xtit = [' ', ' ', ' ', 'Downslope distance (m)']

; define y margins
ymar0 = [2, 4, 6]
ymar1 = [2, 0,-2]
;ymar0 = [2, 4, 6, 8]
;ymar1 = [2, 0,-2, -4]

; define experiments
cExp = ['-exp1','-exp2']

iCol = [80,250]

; define prefix and suffix
fPref = 'syntheticHillslope'
fSuff = '_spinup_testSumma'

; loop through output times
for jTime=0,n_elements(outTime)-1 do begin

 ; define desired time index
 iTime = outTime[jTime]

 ; make a base plot
 plot, indgen(5), xrange=[0,50], yrange=[0,0.5], xstyle=1, ystyle=1, $
  xtitle=xTit[jTime], ytitle=yTit[jTime], ymargin=[ymar0[jTime],yMar1[jTime]], $
  /nodata

 ; plot a legend
 yoff = 0.01
 plots, indgen(6)+30, replicate(0.45,10), psym = sym(6)
 plots, [30,35], [0.4,0.4]
 xyouts, 36, 0.45-yoff, 'SUMMA', charsize=1.5
 xyouts, 36, 0.4-yoff, 'Kinematic', charsize=1.5

 ; plot the time
 xyouts, 2, 0.45, strtrim(iTime,2) + ' hours', charsize=1.5

 ; loop through experiments
 for iExp=0,n_elements(cExp)-1 do begin

  ; define parameters for each experiment
  if(cExp[iExp] eq '-exp1')then begin
   ; experiment 1
   TOPMODEL_exp = 1.0d   ; -
   surfHydCond  = 0.3d   ; m h-1
  endif else begin
   ; experiment 2
   TOPMODEL_exp = 3.0d   ; -
   surfHydCond  = 3.0d   ; m h-1
  endelse
 
  ; define filename
  filenm = fPref + cExp[iExp] + fSuff + '.nc'

  ; open file
  nc_file = ncdf_open(fpath+filenm, /nowrite)

  ; extract the time vector
  ivar_id = ncdf_varid(nc_file,'time')
  ncdf_varget, nc_file, ivar_id, atime

  ; define number of time elements
  ntime = n_elements(atime)-1

  ; loop through HRUs
  for iHRU=0,nHRU-1 do begin

   ; get the number of soil layers
   ivar_id = ncdf_varid(nc_file,'nSoil')
   ncdf_varget, nc_file, ivar_id, nSoil, offset=[iHRU,0], count=[1,ntime]

   ; get the field capacity
   ivar_id = ncdf_varid(nc_file,'fieldCapacity')
   ncdf_varget, nc_file, ivar_id, fieldCapacity, offset=iHRU, count=1

   ; get the start index for mLayers
   ivar_id = ncdf_varid(nc_file,'midSoilStartIndex')
   ncdf_varget, nc_file, ivar_id, midSoilStartIndex, offset=[iHRU,0], count=[1,ntime]

   ; get the depth of each layer
   ivar_id = ncdf_varid(nc_file,'mLayerDepth')
   ncdf_varget, nc_file, ivar_id, mLayerDepth, offset=[iHRU,midSoilStartIndex[0,itime-1]-1], count=[1,nSoil[0,itime-1]]

   ; get the desired variable for all layers
   ivar_id = ncdf_varid(nc_file,'mLayerVolFracLiq')
   ncdf_varget, nc_file, ivar_id, mLayerVolFracLiq, offset=[iHRU,midSoilStartIndex[0,itime-1]-1], count=[1,nSoil[0,itime-1]]

   ; get storage (m)
   totalStorage[iHRU] = total(mLayerDepth[0,*]*(mLayerVolFracLiq[0,*]-fieldCapacity))

  endfor  ; looping through HRUs

  ; define the length of the hillslope at steady state (m)
  t0 = surfHydCond*soilDepth/TOPMODEL_exp         ; m2 h-1
  x1 = t0*tan_slope/qRain                         ; m
  x2 = double(itime)*qRain/(porosity*soilDepth)   ; -
  xL = min([x1*x2^TOPMODEL_exp, totalLength])     ; m

  ; define the steady-state water table
  hS = soilDepth*(qRain*xCoord/(T0*tan_slope))^(1.d/TOPMODEL_exp)

  ; define the water table at length xL
  hL = soilDepth*(qRain*xL/(T0*tan_slope))^(1.d/TOPMODEL_exp)

  ; define the water table from the kinematic predictions
  ix = where(hS lt hL, complement=jx)
  hK = dblarr(nCoord)
  hK[ix] = hS[ix]
  hk[jx] = hL

  ; plot the kinematic predictions
  oplot, xCoord, hK*porosity, color=iCol[iExp]

  ; overplot the model predictions
  oplot, xSlope, totalStorage, color=iCol[iExp], psym=sym(6)

  ; close NetCDF file
  ncdf_close, nc_file

 endfor  ; looping through experiments
endfor  ; looping through output times


; *************************************************************************************************************
; * PART 2: PLOT THE HYDROGRAPHS
; *************************************************************************************************************

; define number of time steps 
ntime = 1000

; define the time coordinate
xTime = dindgen(ntime)+1.d

; define forcing
qrain = dblarr(ntime)
qRain[  0:549] = 0.002  ; m h-1
qRain[550:999] = 0.d

; define the hydrograph from the kinematic predictions
qKin = dblarr(ntime)

; define y margins
ymar0 = [-4.5,-9.5]
ymar1 = [ 2.0, 7.0]

; define x-title
xtit = [' ', 'Time (hours)']

; loop through experiments
for iExp=0,n_elements(cExp)-1 do begin

 ; make a base plot for runoff
 plot, indgen(5), xrange=[0,1000], yrange=[0,2.5], xstyle=1, ystyle=1, $
  xtitle=xtit[iExp], ytitle='Runoff (mm/hour)', $
  ymargin=[ymar0[iExp],ymar1[iExp]], /nodata

 ; plot a legend
 yoff = 0.03
 plots, [600,750], [2.25,2.25], color=iCol[iExp], thick=3
 plots, [600,750], [2.00,2.00]
 xyouts, 760, 2.25-yoff, 'SUMMA', charsize=1.5
 xyouts, 760, 2.00-yoff, 'Kinematic', charsize=1.5

 ; ***
 ; get summa simulations...

 ; define filename
 filenm = fPref + cExp[iExp] + fSuff + '.nc'

 ; open file
 nc_file = ncdf_open(fpath+filenm, /nowrite)

 ; extract the time vector
 ivar_id = ncdf_varid(nc_file,'time')
 ncdf_varget, nc_file, ivar_id, atime

 ;  get the basin column outflow (m3 s-1)
 ivar_id = ncdf_varid(nc_file,'basin__ColumnOutflow')
 ncdf_varget, nc_file, ivar_id, basin__ColumnOutflow

 ; get the total basin area (m2)
 ivar_id = ncdf_varid(nc_file,'basin__totalArea')
 ncdf_varget, nc_file, ivar_id, basin__totalArea

 ; close NetCDF file
 ncdf_close, nc_file

 ; define the model simulations (m3/s --> mm/hour)
 qSim = 1000.d*(3600.d*basin__ColumnOutflow/basin__totalArea)

 ; plot the model simulations
 oplot, xTime[0:nTime-1], qSim[0:ntime-1], linestyle=0, thick=3, color=iCol[iExp]
 
 ; ***
 ; get kinematic predictions...

 ; define parameters for each experiment
 if(cExp[iExp] eq '-exp1')then begin
  ; experiment 1
  TOPMODEL_exp = 1.0d   ; -
  surfHydCond  = 0.3d   ; m h-1
 endif else begin
  ; experiment 2
  TOPMODEL_exp = 3.0d   ; -
  surfHydCond  = 3.0d   ; m h-1
 endelse

 ; define terms
 ; NOTE: state equation:
 ;  porosity*(dz/dt) = dT/dx   (T = transmissivity (m2 h-1) = tan_slope*T0*(z/soilDepth)^TOPMODEL_exp)
 ; can be written written (using the chain rule) as
 ;  porosity*(dz/dt) = (dz/dx) * (dT/dz)
 ; so we need the derivative dT/dz
 ; NOW: the celerity is then
 ;  dx/dt = (1/porosity)*dT/dz
 ; which can be couched in difference form (Beven, 1981) as
 ;  deltaX = deltaTime * dT/dz
 ; and setting deltaX to be the distance from the upstream boundary
 ; the eqn can be used to estimate that the time that the water table height
 ; reaches the downstream boundary 

 ; define the steady-state water table
 t0 = surfHydCond*soilDepth/TOPMODEL_exp      ; m2 h-1
 hS = soilDepth*(qRain[0]*xCoord/(T0*tan_slope))^(1.d/TOPMODEL_exp)

 ; rising limb: calculate q based on water table depth at the downslope end
 for itime=0,ntime-1 do begin

  if(qRain[itime] gt 0.d)then begin
   ; get the water table depth at the downslope end
   tc = xtime[itime]  ; time
   t0 = surfHydCond*soilDepth/TOPMODEL_exp      ; m2 h-1
   x1 = t0*tan_slope/qRain[itime]               ; m
   x2 = tc*qRain[itime]/(porosity*soilDepth)    ; -
   xL = min([x1*x2^TOPMODEL_exp, totalLength])  ; m
   hL = soilDepth*(qRain[itime]*xL/(T0*tan_slope))^(1.d/TOPMODEL_exp)
   ; get the streamflow (m h-1)
   qKin[itime] = (1.d/domainArea)*hillwidth*tan_slope*t0*(hL/soilDepth)^TOPMODEL_exp

  ; falling limb 
  endif else begin

   qKin[itime] = -9999.d

  endelse

 endfor  ; (looping through time)

 ; falling limb: calculate q based on the time that it takes for the water table to reach the downstream end

 ; compute celerity of the water table for all positions on the hillslope; dxdt = f(hS)
 t0 = surfHydCond*soilDepth/TOPMODEL_exp      ; maximum transmissivity m2 h-1
 dtdz = (tan_slope*t0/soilDepth)*TOPMODEL_exp*(hS/soilDepth)^(TOPMODEL_exp - 1.d)  ; derivative in transmissivity w.r.t water table height (m h-1)
 dxdt = dtdz/porosity  ; celerity (m h-1)
 ; define the time that a given water table height will reach the downstream boundary
 td   = 550.d                 ; time when drying begins
 delX = totalLength - xCoord  ; distance from the downstream end
 tExit = td + delX/dxdt        ; time that a water table height will reach the downstream end (h-1
 ; compute the flow associated with the water table height
 qFall = (1.d/domainArea)*hillwidth*tan_slope*t0*(hS/soilDepth)^TOPMODEL_exp

 ; plot the kinematic predictions
 oplot, xTime, qKin*1000.d, color=0, min_value = 0.d  ; rising limb
 oplot, tExit, qFall*1000.d, color=0  ; falling limb

endfor ; looping through experiments

; make a figure
write_png, 'zFigures/wigmosta1999_testSUMMA.png', tvrd(true=1)

stop
end
