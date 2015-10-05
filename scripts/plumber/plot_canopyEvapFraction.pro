pro plot_canopyEvapFraction

; define plotting parameters
window, 0, xs=2050, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3.5
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,1]

; define the path to the plumber data
plumber_path = '/d1/mclark/PLUMBER_data/'

; define the path to the model output
model_path = plumber_path + 'model_output/'

; define the path to the site data
site_path =  plumber_path + 'site_data/'

; define the model names
model_names = [$
               'SUMMA.1.0.exp.02.022', $
               'SUMMA.1.0.exp.02.017', $
               'SUMMA.1.0.exp.02.021', $
               'SUMMA.1.0.exp.02.024', $
               'SUMMA.1.0.exp.02.025', $
               'SUMMA.1.0.exp.02.028', $
               'SUMMA.1.0.exp.02.029', $
               ' ']
nModels = n_elements(model_names)-1  ; -1 to take care of the blank at the end

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
nSites  = n_elements(site_names)

; define desired parameters
cDesireParam = ['refInterceptCapRain', $
                'canopyWettingFactor', $
                'canopyWettingExp']
nParam = n_elements(cDesireParam)

; *****
; * GET DATA...
; *************

; skip if got data already
goto, got_data

; define an array for the model data
canLatHeat = dblarr(nModels,nSites)
totLatHeat = dblarr(nModels,nSites)

; define vectors for the observations
Precip = dblarr(nSites)
QleObs = dblarr(nSites)

; define vectors for the parameters
xParam = dblarr(nParam,nModels)

; loop through sites
for iSite=0,nSites-1 do begin

 ; get the precip data
 nc_file = ncdf_open(site_path + 'met/' + site_names[iSite] + 'Fluxnet.1.4_met.nc')
  ivar_id = ncdf_varid(nc_file,'Rainf')
  ncdf_varget, nc_file, ivar_id, RainfVector
  Precip[iSite] = mean(RainfVector)
 ncdf_close, nc_file

 ; get the observed latent heat flux
 nc_file = ncdf_open(site_path + 'flux/' + site_names[iSite] + 'Fluxnet.1.4_flux.nc')
  ivar_id = ncdf_varid(nc_file,'Qle')
  ncdf_varget, nc_file, ivar_id, QleVector
  QleObs[iSite] = mean(QleVector)
 ncdf_close, nc_file

 ; loop thru models
 for iModel=0,nModels-1 do begin

  ; define filename
  filename = model_path + model_names[iModel] + '/' + model_names[iModel] + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'

  ; open file for reading
  nc_file = ncdf_open(filename, /nowrite)

   ; get the parameters
   if(iSite eq 0)then begin
    for iParam=0,nParam-1 do begin
     ivar_id = ncdf_varid(nc_file,cDesireParam[iParam])
     ncdf_varget, nc_file, ivar_id, xTemp
     xParam[iParam,iModel] = reform(xTemp)
    endfor
   endif

   ; read in canopy evap
   ivar_id = ncdf_varid(nc_file,'scalarLatHeatCanopyEvap')
   ncdf_varget, nc_file, ivar_id, xTemp
   canLatHeat[iModel,iSite] = -mean(xTemp)  ; NOTE: reverse sign

   ; read in latent heat
   ivar_id = ncdf_varid(nc_file,'Qle')
   ncdf_varget, nc_file, ivar_id, xTemp
   totLatHeat[iModel,iSite] = -mean(xtemp)  ; NOTE: reverse sign

  ; close netcdf file
  ncdf_close, nc_file

  ; print progress
  print, filename, n_elements(xTemp), n_elements(RainfVector)

 endfor  ; looping thru realizations
 
 ; print break
 print, '**'

endfor  ; looping through sites

; save statistics
save, xparam, canLatHeat, totLatHeat, QleObs, Precip, filename='xIDLsave/canEvapFraction.sav'

; label if have the statistics already
got_data:

; restore the statistics
restore, 'xIDLsave/canEvapFraction.sav'

; *****
; * MAKE SCATTER PLOTS...
; ***********************
; define the latent heat of vaporization
LHvap = 2501000.d  ; J kg-1

; define x-margins
xmar1 = [8,8]
xmar2 = [1,1]

; define y-margin
ymar = [4,1]

; define colors
ixColor=[80,160,250,80,160,210,250]

; define symbols
ixSymbol = [6,6,6,1,1,1,1]

; get the ratio of evaporation to precipitation
eOverP_obs = reform((QleObs[*]/LHvap)/Precip[*])

; define arrays for the model simulations
cOverT_mod = fltarr(nModels,nSites)
eOverP_mod = fltarr(nModels,nSites)

; loop through models
for iModel=0,nModels-1 do begin

 ; get the ratio of evaporation to precipitation
 eOverP_mod[iModel,*] = reform((totLatHeat[iModel,*]/LHvap)/Precip[*])

 ; get the ratio of canopy evap to total evap
 cOverT_mod[iModel,*] = reform(canLatHeat[iModel,*]) / reform(totLatHeat[iModel,*]) 

endfor  ; looping through models

; define the base simulation
ixBase = 1

; loop through plots
for jplot=0,1 do begin

 ; define the margins
 xmar=[xmar1[jplot],xmar2[jplot]]

 ; define xmax
 if(jplot eq 0)then xmax=0.3 else xmax = 1.15

 ; define ymax
 if(jplot eq 0)then ymax=0.525 else ymax = 1.15

 ; define the x and y title
 if(jplot eq 0)then xtitle='Canopy evaporation fraction (parset 2)' else xtitle='Observed precipitation fraction' 
 if(jplot eq 0)then ytitle='Canopy evaporation fraction (other parsets)' else ytitle='SUMMA precipitation fraction' 

 ; make a base plot
 plot, indgen(5), xrange=[0,xmax], yrange=[0,ymax], xstyle=1, ystyle=1, $
  xtitle=xtitle, ytitle=ytitle, xmargin=xmar, ymargin=ymar, /nodata
 plots, [0,xmax], [0,xmax]

 ; loop through models
 for iModel=0,nModels-1 do begin
  if(jplot eq 0)then begin
   if(iModel ne ixBase)then oplot, cOverT_mod[ixBase,*], cOverT_mod[iModel,*], psym=sym(ixSymbol[iModel]), color=ixColor[iModel], symsize=5
  endif else begin
   oplot, eOverP_obs[*], eOverP_mod[iModel,*], psym=sym(ixSymbol[iModel]), color=ixColor[iModel], symsize=5
  endelse
 endfor  ; looping thru models

 ; plot a legend
 ytop = 0.30*ymax
 yInc = 0.04*ymax
 iLegend=0
 for jLegend=0,6 do begin
  if(jplot eq 0 and jLegend eq 1)then continue
  ypos = ytop - iLegend*yInc
  legendText = 'Parset ' + strtrim(jLegend+1,2)
  plots, 0.75*xmax, ypos, psym=sym(ixSymbol[jLegend]), color=ixColor[jLegend], symsize=5
  xyouts, 0.77*xmax, ypos-0.015*ymax, legendText
  iLegend = iLegend+1
 endfor

endfor ; looping through plot types

; write figure
write_png, 'figures/summaCanEvapFraction.png', tvrd(true=1)

stop
end
