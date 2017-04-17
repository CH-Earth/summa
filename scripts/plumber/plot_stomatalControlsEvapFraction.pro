pro plot_stomatalControlsEvapFraction

; define plotting parameters
window, 0, xs=1350, ys=1200, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2.5
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,2]

; define the path to the plumber data
plumber_path = '/d1/mclark/PLUMBER_data/'

; define the path to the model output
model_path = plumber_path + 'model_output/'

; define the path to the site data
site_path =  plumber_path + 'site_data/'

; define the model names
model_names = [$
               'SUMMA.1.0.exp.02.013', $
               'SUMMA.1.0.exp.02.014', $
               'SUMMA.1.0.exp.02.027', $
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
ixColor=[80,160,250]

; define symbols
ixSymTyp = [5,1,12]

; define symbol size
ixSymsize = [5,4,3]

; define legend
legendText = ['Ball-Berry (h1)', $
              'Ball-Berry (h2)', $
              'Jarvis']

; define axis
axisName = ['[Ball-Berry: linear humidity function]', $
            '[Ball-Berry: hyperbolic humidity function]', $
            '[Jarvis]']

; get the ratio of evaporation to precipitation
eOverP_obs = reform((QleObs[*]/LHvap)/Precip[*])

; define arrays for the model simulations
eOverP_mod = fltarr(nModels,nSites)

; get the ratio of evaporation to precipitation
for iModel=0,nModels-1 do begin
 eOverP_mod[iModel,*] = reform((totLatHeat[iModel,*]/LHvap)/Precip[*])
endfor  ; looping through models

; define the x and y title
xtitle='Observed precipitation fraction' 
ytitle='SUMMA precipitation fraction' 

; define xmax and ymax
xmax=1.5
ymax=1.5

; plot model vs model
for jModel=1,nModels-1 do begin
 for imodel=0,jModel-1 do begin
  ; define axis titles
  xtitle = 'Evaporative fraction!C' + axisName[iModel]
  ytitle = 'Evaporative fraction!C' + axisName[jModel]
  ; make a base plot
  plot, indgen(5), xrange=[0,xmax], yrange=[0,ymax], xstyle=1, ystyle=1, $
   xtitle=xtitle, ytitle=ytitle, xmargin=xmar, ymargin=ymar, /nodata
  plots, [0,xmax], [0,xmax]
  ; plot scatter
  for iSite=0,nSites-1 do begin
   plots, eOverP_mod[iModel,iSite], eOverP_mod[jModel,iSite], psym=sym(1), color=250, symsize=4
   plots, eOverP_mod[iModel,iSite], eOverP_mod[jModel,iSite], psym=sym(6), color=0, symsize=4
  endfor
 endfor
endfor  ; looping thru models

; make a base plot
plot, indgen(5), xrange=[0,xmax], yrange=[0,ymax], xstyle=1, ystyle=1, $
 xtitle='Observed evaporative fraction', ytitle='SUMMA evaporative fraction', $
 xmargin=xmar, ymargin=ymar, /nodata
plots, [0,xmax], [0,xmax]

; plot obs vs model
for iModel=0,nModels-1 do begin
 oplot, eOverP_obs[*], eOverP_mod[iModel,*], psym=sym(ixSymTyp[iModel]), color=ixColor[iModel], symsize=ixSymSize[iModel], thick=2
endfor

; plot a legend
ytop = 0.24*ymax
yInc = 0.075*ymax
for jLegend=0,n_elements(legendText)-1 do begin
 ypos = ytop - jLegend*yInc
 plots, 0.53*xmax, ypos, psym=sym(ixSymTyp[jLegend]), color=ixColor[jLegend], symsize=ixSymSize[jLegend], thick=2
 xyouts, 0.57*xmax, ypos-0.015*ymax, legendText[jLegend]
endfor


; write figure
write_png, 'figures/stomatalControlsEvapFraction.png', tvrd(true=1)

stop
end
