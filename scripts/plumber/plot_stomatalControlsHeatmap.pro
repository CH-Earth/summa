pro plot_stomatalControlsHeatmap

; define plotting parameters
window, 0, xs=1300, ys=1200, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=1.5
!P.COLOR=1
erase, color=255
!P.MULTI=[0,2,2,0,0]

; refine colors
tvlct, r, g, b, /get
r[0] = 255
g[0] = 255
b[0] = 255
tvlct, r, g, b

; define file path
fPath = '/d1/mclark/PLUMBER_data/model_output/'

; define the site names
site_names = ['Amplero',     $
              'Blodgett',    $
              'Bugac',       $
              'ElSaler2',    $
              'ElSaler',     $
              'Espirra',     $
              'FortPeck',    $
              'Harvard',     $
              'Hesse',       $
              'Howard',      $
              'Howlandm',    $
              'Hyytiala',    $
              'Kruger',      $
              'Loobos',      $
              'Merbleue',    $
              'Mopane',      $
              'Palang',      $
              'Sylvania',    $
              'Tumba',       $
              'UniMich'      ]

; define the number of sites
nSites  = n_elements(site_names)

;ix        vegeParTbl        soilStress        stomResist        bbTempFunc        bbHumdFunc        bbElecFunc        bbAssimFnc        bbCanIntg8        rootProfil
;-- ----------------- ----------------- ----------------- ----------------- ----------------- ----------------- ----------------- ----------------- -----------------
; 0              USGS          NoahType         BallBerry           q10Func  humidLeafSurface            linear           minFunc   constantScaling         doubleExp
; 1              USGS          NoahType     BallBerryFlex           q10Func  humidLeafSurface            linear           minFunc   constantScaling         doubleExp
; 2              USGS          NoahType     BallBerryFlex           q10Func  humidLeafSurface            linear           minFunc   constantScaling         doubleExp
; 3              USGS          NoahType     BallBerryFlex           q10Func  scaledHyperbolic            linear           minFunc   constantScaling         doubleExp
; 4              USGS          NoahType     BallBerryFlex         Arrhenius  humidLeafSurface            linear           minFunc   constantScaling         doubleExp
; 5              USGS          NoahType     BallBerryFlex         Arrhenius  scaledHyperbolic            linear           minFunc   constantScaling         doubleExp
; 6              USGS          NoahType     BallBerryFlex         Arrhenius  humidLeafSurface            linear      colimitation   constantScaling         doubleExp
; 7              USGS          NoahType     BallBerryFlex         Arrhenius  scaledHyperbolic            linear      colimitation   constantScaling         doubleExp
; 8              USGS          NoahType     BallBerryFlex         Arrhenius  humidLeafSurface     quadraticJmax      colimitation   constantScaling         doubleExp
; 9              USGS          NoahType     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation   constantScaling         doubleExp
;10              USGS          NoahType     BallBerryFlex         Arrhenius  humidLeafSurface     quadraticJmax      colimitation        laiScaling         doubleExp
;11              USGS          NoahType     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp
;12              USGS          NoahType            Jarvis         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp
;13              USGS          NoahType     BallBerryFlex         Arrhenius  humidLeafSurface     quadraticJmax      colimitation        laiScaling          powerLaw
;14              USGS          NoahType     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling          powerLaw
;15              USGS          NoahType     BallBerryFlex         Arrhenius  humidLeafSurface     quadraticJmax      colimitation        laiScaling          powerLaw
;16              USGS          NoahType     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling          powerLaw
;17              USGS          CLM_Type     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp
;18              USGS          CLM_Type     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation   constantScaling         doubleExp
;19              USGS          SiB_Type     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp
;21              USGS          CLM_Type     BallBerryFlex         Arrhenius  humidLeafSurface     quadraticJmax      colimitation        laiScaling         doubleExp
;22              USGS          CLM_Type     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp
;25              USGS          CLM_Type     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp
;26              USGS          CLM_Type     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp
;30      plumberCABLE          CLM_Type     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp
;31   plumberCHTESSEL          CLM_Type     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp
;32      plumberSUMMA          CLM_Type     BallBerryFlex         Arrhenius  scaledHyperbolic     quadraticJmax      colimitation        laiScaling         doubleExp

; define simulation pairs
simPairs = [ ['002', '004'], $      ; q10 vs Arrhenius (humidLeafSurface, linear, minFunc, constantScaling, doubleExp)
             ['006', '008'], $      ; linear vs quadraticJmax (Arrhenius, humidLeafSurface, colimitation, constantScaling, doubleExp)
             ['010', '011'], $      ; humidLeafSurface vs scaledHyperbolic (Arrhenius, quadraticJmax, laiScaling, doubleExp)
             ['008', '010']  ]      ; constant vs LAI scaling (Arrhenius, humidLeafSurface, quadraticJmax, doubleExp)

; define the title
titleVec = ['Temperature function', $
            'PAR function',      $
            'Humidity function',    $
            'Canopy integration']

; define x title
xtitleVec = ['Stomatal conductance (m s!e-1!n)!C' + '['+titleVec[0] + ' = ' + 'q10'     + ']', $
             'Stomatal conductance (m s!e-1!n)!C' + '['+titleVec[1] + ' = ' + 'Linear'  + ']', $
             'Stomatal conductance (m s!e-1!n)!C' + '['+titleVec[2] + ' = ' + 'Linear'  + ']', $
             'Stomatal conductance (m s!e-1!n)!C' + '['+titleVec[3] + ' = ' + 'Constant'+ ']']

; define y title
ytitleVec = ['Stomatal conductance (m s!e-1!n)!C' + '['+titleVec[0] + ' = ' + 'Arrhenious'        + ']', $
             'Stomatal conductance (m s!e-1!n)!C' + '['+titleVec[1] + ' = ' + 'Quadratic Jmax'    + ']', $
             'Stomatal conductance (m s!e-1!n)!C' + '['+titleVec[2] + ' = ' + 'Hyperbolic'        + ']', $
             'Stomatal conductance (m s!e-1!n)!C' + '['+titleVec[3] + ' = ' + 'LAI-dependent'     + ']'  ]

; get size
iSize = size(simPairs, /dimensions)
nCompare = iSize[1]

; define the heat map
dx       = 0.2
dy       = 0.2
nx       = 100
ny       = 100
xHeatMap = dx*dindgen(nx)
yHeatMap = dy*dindgen(ny)
zHeatMap = dblarr(nx,ny)

; define the x margin
xmar1 = [12, 6,12, 6]
xmar2 = [10,16,10,16]

; loop through comparisons
for iCompare=0,nCompare-1 do begin

 ; define models
 models = simPairs[*,iCompare]

 ; define the heat map file 
 heatMapSaveFile = 'xIDLsave/stomatalControlHeatmap_'+models[0]+'_'+models[1]+'.sav'

 ; define titles
 xtitle=xtitleVec[iCompare]
 ytitle=ytitleVec[iCompare]

 ; define margins
 xmar=[xmar1[iCompare],xmar2[iCompare]]
 ymar=[8,2]

 ; make a base plot
 plot, xHeatMap, yHeatMap, xrange=[xHeatMap[0],xHeatMap[nx-1]+dx], yrange=[yHeatMap[0],yHeatMap[ny-1]+dy], $
  xstyle=9, ystyle=9, xtitle=xtitle, ytitle=ytitle, xmargin=xmar, ymargin=ymar, $
  xcharsize=1.5, ycharsize=1.5, xticklen=(-0.02), yticklen=(-0.02), /nodata

 ; check if the file exists
 if(file_test(heatMapSaveFile) eq 0) then begin

  ; initialize the heat map
  zHeatMap[*,*] = 0.d

  ; loop through sites
  for iSite=0,nSites-1 do begin

   ; print progress
   print, 'processing site ' + site_names[isite]

   ; read in data
   for ifile=0,1 do begin

    ; define model
    cModel = 'SUMMA.1.0.exp.02.' + models[ifile]

    ; define file
    filename = fPath + cModel + '/' + cModel + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'

    ; open file
    nc_file = ncdf_open(filename, /nowrite)

    ; read total stomatal conductance
    ivar_id = ncdf_varid(nc_file,'stomatalConductance')
    ncdf_varget, nc_file, ivar_id, totalStomatalConductance

    ; save data
    if(ifile eq 0)then begin
     total_g01    = reform(totalStomatalConductance)*1000.d
    endif else begin
     total_g02    = reform(totalStomatalConductance)*1000.d
    endelse

   endfor  ; looping through files

   ; get the heat map
   for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
     iMatch = where(total_g01 gt xHeatMap[ix] and total_g01 le xHeatMap[ix]+dx and $
                    total_g02 gt yHeatMap[iy] and total_g02 le yHeatMap[iy]+dy, nMatch)
     zHeatMap[ix,iy] = zHeatMap[ix,iy] + double(nMatch)
    endfor
   endfor

   ; plot the heat map
   for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
     xx = [xHeatMap[ix], xHeatMap[ix]+dx, xHeatMap[ix]+dx, xHeatMap[ix]]
     yy = [yHeatMap[iy], yHeatMap[iy], yHeatMap[iy]+dy, yHeatMap[iy]+dy]
     if(zHeatMap[ix,iy] gt 0.)then zz = min([alog10(zHeatMap[ix,iy])*75.d, 250.d]) else zz = 0.d
     polyfill, xx, yy, color=zz
    endfor
   endfor

  endfor  ; looping through sites

  ; save the heat map
  save, zHeatMap, filename=heatMapSaveFile

 ; heat map save file exists
 endif else begin

  ; restore the statistics
  restore, heatMapSaveFile 

  ; plot the heat map
  for ix=0,nx-1 do begin
   for iy=0,ny-1 do begin
    xx = [xHeatMap[ix], xHeatMap[ix]+dx, xHeatMap[ix]+dx, xHeatMap[ix]]
    yy = [yHeatMap[iy], yHeatMap[iy], yHeatMap[iy]+dy, yHeatMap[iy]+dy]
    if(zHeatMap[ix,iy] gt 0.)then zz = min([alog10(zHeatMap[ix,iy])*75.d, 250.d]) else zz = 0.d
    polyfill, xx, yy, color=zz
   endfor
  endfor

  ; re-draw axes
  plots, [0, xHeatMap[nx-1]+dx], [0, 0]
  plots, [0, xHeatMap[nx-1]+dx], [yHeatMap[ny-1]+dy, yHeatMap[ny-1]+dy]
  plots, [xHeatMap[nx-1]+dx, xHeatMap[nx-1]+dx], [0,yHeatMap[ny-1]+dy]
  plots, [0, 0], [0,yHeatMap[ny-1]+dy]

  ; plot 1:1 line
  plots, [0,xHeatMap[nx-1]+dx], [0,yHeatMap[ny-1]+dy]
 
  ; plot the title
  xyouts, 0.05*(xHeatMap[nx-1]+dx), 0.90*(yHeatMap[ny-1]+dy), titleVec[iCompare], charsize=3

 endelse  ; if heat map save file exists

endfor  ; looping through comparisons

; plot a colorbar
ypos1 = 0.10
ypos2 = 0.97
xColors = indgen(251)
xLabels = xColors/75.d
colorbar, 0.91, 0.94, ypos1, ypos2, xLabels, xColors, every=15, charsize=1.5, /nobox, /norm
xyouts, 0.985, 0.5*(ypos1 + ypos2), 'log!i10!n count', alignment=0.5, orientation=90, charsize=3, /normal

; write figure
write_png, 'figures/stomatalControls.png', tvrd(true=1)

stop
end
