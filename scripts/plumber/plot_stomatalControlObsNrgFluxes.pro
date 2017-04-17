pro plot_stomatalControlObsNrgFluxes

; define plotting parameters
window, 0, xs=2050, ys=1200, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3
!P.COLOR=1
erase, color=255
!P.MULTI=[0,3,2,0,0]

; refine colors
tvlct, r, g, b, /get
r[0] = 255
g[0] = 255
b[0] = 255
tvlct, r, g, b

; define the path to the plumber data
plumber_path = '/d1/mclark/PLUMBER_data/'

; define the path to the model output
model_path = plumber_path + 'model_output/'

; define the path to the site data
site_path =  plumber_path + 'site_data/flux/'

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
simPairs = [ ['obs', '014'], $      ; Ball-Berry vs obs
             ['obs', '027'], $      ; Jarvis vs obs
             ['014', '027']  ]      ; Ball-Berry vs Jarvis

; define variables
varDesire = ['Qh','Qle']

varTitle = ['Sensible heat (W m!e-2!n)','Latent heat (W m!e-2!n)']

; define x title
xtitleVec = ['Observations', $
             'Observations', $
             'SUMMA [Ball-Berry]']

; define y title
ytitleVec = ['SUMMA [Ball-Berry]', $
             'SUMMA [Jarvis]', $
             'SUMMA [Jarvis]']

; get size
iSize = size(simPairs, /dimensions)
nCompare = iSize[1]

; define the heat map
dx       = 7.5
dy       = 7.5
nx       = 100
ny       = 100
xHeatMap = dx*dindgen(nx) - 0.2*dx*double(nx)
yHeatMap = dy*dindgen(ny) - 0.2*dy*double(ny)
zHeatMap = dblarr(nx,ny)

; define the x margin
xmar1 = [16, 8, 0]
xmar2 = [10,18,26]

; loop through variables
for iVar=0,n_elements(varDesire)-1 do begin 

 ; loop through comparisons
 for iCompare=0,nCompare-1 do begin

  ; define models
  models = simPairs[*,iCompare]

  ; define the heat map file 
  heatMapSaveFile = 'xIDLsave/controlsNrgFluxesHeatmap_'+models[0]+'_'+models[1]+'_'+varDesire[iVar]+'.sav'
  print, 'heatMapSaveFile = ', heatMapSaveFile

  ; define titles
  xtitle=varTitle[iVar] + '!C' + xtitleVec[iCompare]
  ytitle=varTitle[iVar] + '!C' + ytitleVec[iCompare]

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

     ; define file name
     if(models[ifile] eq 'obs')then begin
      filename = site_path + site_names[iSite] + 'Fluxnet.1.4_flux.nc'
     endif else begin
      cModel   = 'SUMMA.1.0.exp.02.' + models[ifile]
      filename = model_path + cModel + '/' + cModel + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'
     endelse
     print, '   - filename = ', filename

     ; open file
     nc_file = ncdf_open(filename, /nowrite)

     ; read desired variable
     ivar_id = ncdf_varid(nc_file,varDesire[iVar])
     ncdf_varget, nc_file, ivar_id, tmpdat

     ; save data
     if(ifile eq 0)then begin
      if(models[ifile] eq 'obs')then dat01=reform(tmpdat) else dat01= -reform(tmpdat)
     endif else begin
      if(models[ifile] eq 'obs')then dat02=reform(tmpdat) else dat02= -reform(tmpdat)
     endelse

    endfor  ; looping through files

    ; get the heat map
    for ix=0,nx-1 do begin
     for iy=0,ny-1 do begin
      iMatch = where(dat01 gt xHeatMap[ix] and dat01 le xHeatMap[ix]+dx and $
                     dat02 gt yHeatMap[iy] and dat02 le yHeatMap[iy]+dy, nMatch)
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

   ; define points
   x0 = xHeatMap[0]
   x1 = xHeatMap[nx-1]+dx
   y0 = yHeatMap[0]
   y1 = yHeatMap[ny-1]+dy

   ; re-draw axes
   plots, [x0,x1], [y0,y0]  ; bottom
   plots, [x1,x1], [y0,y1]  ; right
   plots, [x0,x1], [y1,y1]  ; top
   plots, [x0,x0], [y0,y1]  ; left

   ; plot 1:1 line
   plots, [x0,x1], [y0,y1]
 
  endelse  ; if heat map save file exists

 endfor  ; looping through comparisons
endfor  ; looping through the variables

; plot a colorbar
ypos1 = 0.10
ypos2 = 0.97
xColors = indgen(251)
xLabels = xColors/75.d
colorbar, 0.91, 0.94, ypos1, ypos2, xLabels, xColors, every=15, charsize=2.5, /nobox, /norm
xyouts, 0.985, 0.5*(ypos1 + ypos2), 'log!i10!n count', alignment=0.5, orientation=90, charsize=4, /normal

; write figure
write_png, 'figures/stomatalControlsObsNrgFluxes.png', tvrd(true=1)

stop
end
