pro plot_budyko

; define plotting parameters
window, 0, xs=1500, ys=500, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
;!P.MULTI=[0,2,2,0,1]
!P.MULTI=1

; define the path to the plumber data
plumber_path = '/Volumes/d1/mclark/PLUMBER_data/'

; define the path to the model output
model_path = plumber_path + 'model_output/'

; define the path to the site data
site_path =  plumber_path + 'site_data/'

; define the benchmark path
benchmark_path = plumber_path + 'benchmark_data/'

; define the benchmarks
stat_benchmarks = ['_1lin','_2lin','_3km27']
phys_benchmarks = ['Penman_Monteith.1','Manabe_Bucket.2']

; define the model names
model_names = ['CABLE.2.0',                  $
               'CABLE_2.0_SLI.vxh599_r553',  $
               'CHTESSEL',                   $
               'COLASSiB.2.0',               $
               'ISBA_SURFEX_3l.SURFEX7.3',   $
               'ISBA_SURFEX_dif.SURFEX7.3',  $
               'JULES.3.1',                  $
               'JULES3.1_altP',              $
               'Mosaic.1',                   $
               'NOAH.2.7.1',                 $
               'Noah.3.2',                   $
               'NOAH.3.3',                   $
               'ORCHIDEE.trunk_r1401'        ] 

; define model symbols
;   1 : filled circle
;   2 : filled upward triangle
;   3 : filled downward triangle
;   4 : filled diamond
;   5 : filled square
;   6 : open circle
;   7 : open upward triangle
;   8 : open downward triangle
;   9 : open diamond
;  10 : open square
;  11 : plus
;  12 : X
;  13 : star
;  16 : open rightfacing triangle
;  17 : open leftfacing triangle
model_symbols = [6,6,7,8,17,17,16,16,13,12,12,12,10,10]
benchmark_symbols = [1,2,3,4,5]

; define a legend
model_legend = ['CABLE v2',  $
                'CABLE-SLI', $
                'CH-TESSEL', $
                'COLA-SSiB', $
                'ISBA-3L',   $
                'ISBA-dif',  $
                'JULES v3',  $
                'JULES altP',$
                'Mosaic',    $
                'Noah v2.7', $
                'Noah v3.2', $
                'Noah v3.3', $
                'ORCHIDEE']
benchmark_legend = ['Penman','Manabe','1 linear','2 linear','3 k-means 27']

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

; define the number of models and sites
nStat   = n_elements(stat_benchmarks)
nPhys   = n_elements(phys_benchmarks)
nModels = n_elements(model_names)
nSites  = n_elements(site_names)

; define the number of realizations
nReal = nModels + nStat + nPhys + 1

; define named variables to define the type of simulation
ixObs   = 3
ixStat  = 0
ixPhys  = 1
ixModel = 2
ixType  = intarr(nReal)

; define type of simulation
ixType[0:nModels-1] = ixModel
ixType[nModels:nModels+nPhys-1] = ixPhys
ixType[nModels+nPhys:nModels+nPhys+nStat-1] = ixStat
ixType[nReal-1] = ixObs

; get the list of models
cList = [model_names,phys_benchmarks,stat_benchmarks,'obs']

; define arrays
Qh   = dblarr(nReal,nSites)
Qle  = dblarr(nReal,nSites)

Precip = dblarr(nSites)

; *****
; * GET DATA...
; *************

; skip if got data already
goto, got_data

; loop through sites
for iSite=0,nSites-1 do begin

 ; initialize counter for the stats varname
 ixStatCount=0

 ; get the precip data
 nc_file = ncdf_open(site_path + 'met/' + site_names[iSite] + 'Fluxnet.1.4_met.nc')
  ivar_id = ncdf_varid(nc_file,'Rainf')
  ncdf_varget, nc_file, ivar_id, RainfVector
  Precip[iSite] = mean(RainfVector)
 ncdf_close, nc_file

 ; loop thru realizations
 for iReal=0,nReal-1 do begin

  ; define filename
  case ixType[iReal] of
   ixObs:   filename = site_path + 'flux/' + site_names[iSite] + 'Fluxnet.1.4_flux.nc'
   ixStat:  filename = benchmark_path + site_names[iSite] + 'Fluxnet_1.4_PLUMBER_benchmarks.nc'
   ixPhys:  filename = model_path + phys_benchmarks[iReal-nModels] + '/' + phys_benchmarks[iReal-nModels] + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'
   ixModel: filename = model_path + model_names[iReal] + '/' + model_names[iReal] + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'
  endcase

  ; define suffix to the variable name
  if(ixType[iReal] eq ixStat)then begin
   vSuffix = stat_benchmarks[ixStatCount]
   ixStatCount = ixStatCount+1
  endif else begin
   vSuffix = ''
  endelse

  ; open file for reading
  nc_file = ncdf_open(filename, /nowrite)

   ; read in sensible heat
   ivar_id = ncdf_varid(nc_file,'Qh'+vSuffix)
   ncdf_varget, nc_file, ivar_id, QhVector
   Qh[iReal,iSite] = mean(QhVector)

   ; read in latent heat
   ivar_id = ncdf_varid(nc_file,'Qle'+vSuffix)
   ncdf_varget, nc_file, ivar_id, QleVector
   Qle[iReal,iSite] = mean(QleVector)

  ; close netcdf file
  ncdf_close, nc_file

  ; print progress
  print, filename, n_elements(QhVector), n_elements(RainfVector)

  ; fix CHTESSEL
  if(ixType[iReal] eq ixModel)then begin
   if(model_names[iReal] eq 'CHTESSEL')then begin
    Qh[iReal,iSite]  = -Qh[iReal,iSite]
    Qle[iReal,iSite] = -Qle[iReal,iSite]
   endif
  endif

 endfor  ; looping thru realizations
 
 ; print break
 print, '**'

endfor  ; looping through sites

; save statistics
save, Qh, Qle, Precip, filename='xIDLsave/budykoData.sav'

; label if have the statistics already
got_data:

; restore the statistics
restore, 'xIDLsave/budykoData.sav'

; *****
; * MAKE BUDYKO PLOTS...
; **********************

; identify where the model is Penman-Monteith
iMatch = where(cList eq 'Penman_Monteith.1', nMatch)
if(nMatch eq 1)then ixMatch=iMatch[0] else stop, 'difficulty identifying PET'

; define the latent heat of vaporization
LHvap = 2501000.d  ; J kg-1

; define xmax
xmax = 10

; define margins
xmar=[8,2]
ymar=[4,1]

; define the x and y title
xtitle='Dryness index'
ytitle='Evaporative fraction'

; define colors
iColor=[250,120,80,0]

; define symbols
iSymType=[1,1,6,1]

; define symbol size
iSymSize=[2,2,1,2]

; define the Turc-Pike curve
nTrial = 1000
vv = -2.d
xx = xmax * (dindgen(nTrial)+1.d)/double(nTrial)
tp = (1.d + xx^vv)^(1.d/vv)

; make a base plot
plot, indgen(5), xrange=[0,xmax], yrange=[0,3.5], xstyle=1, ystyle=1, $ 
 xtitle=xtitle, ytitle=ytitle, xmargin=xmar, ymargin=ymar, /nodata
plots, [0,   1], [0,1]
plots, [1,xmax], [1,1]
oplot, xx, tp

; loop through data types
for jxType=0,3 do begin

 ; loop through realizations
 for iReal=0,nReal-1 do begin

  ; check if matching the current data type
  if(ixType[iReal] ne jxType)then continue

  ; get the available energy
  nrg = reform(Qle[ixMatch,*])

  ; get the dryness index
  dryness = replicate(-9999.d, nSites)
  ixValid = where(nrg gt 0.d, nValid) & if(nValid eq 0)then stop, 'no valid data'
  dryness[ixValid] = (nrg[ixValid]/LHvap)/Precip[ixValid]  ; W m-2 --> kg m-2 s-1

  ; get the ratio of evaporation to precipitation
  eOverP = reform((Qle[iReal,*]/LHvap)/Precip[*]) 

  ; plot the data
  jxColor = iColor[ixType[iReal]]    ; color
  jxSymTy = iSymType[ixType[iReal]]  ; symbol type
  jxSymSz = iSymSize[ixType[iReal]]  ; symbol size
  oplot, dryness, eOverP, psym=sym(jxSymTy), color=jxColor, symsize=jxSymSz

 endfor  ; looping through realizations
endfor  ; looping through data types

; plot a legend
ytop = 3.25
yInc = 0.25
legendText=['Statistical benchmarks','Physical benchmarks','LSM simulations','Observations']
for iLegend=0,3 do begin
 ypos = ytop - iLegend*yInc
 plots, 0.3, ypos, psym=sym(iSymType[iLegend]), color=iColor[iLegend], symsize=iSymSize[iLegend]
 xyouts, 0.4, ypos-0.05, legendText[iLegend]
endfor

; write figure
write_png, 'figures/plumberBudyko.png', tvrd(true=1)

; *****
; * MAKE SCATTER PLOTS...
; ***********************
window, 2, xs=1500, ys=1000, retain=2
device, decomposed=0
!P.MULTI=[0,2,1]

; define xmax and ymax
xmax = 1.5
ymax = 3.0

; define x-margins
xmar1 = [10,5]
xmar2 = [0,5]

; define y-margin
ymar = [4,1]

; define arrays
eOverP_save = dblarr(nReal,nSites)

; get the Budyko prediction
nrg        = reform(Qle[ixMatch,*])      ; energy
dryness    = (nrg/LHvap)/Precip          ; dryness index, W m-2 --> kg m-2 s-1
budykoPred = (1.d + dryness^vv)^(1.d/vv) ; Budyko prediction

; get the ratio of evaporation to precipitation
eOverP_obs = reform((Qle[nReal-1,*]/LHvap)/Precip[*])

; loop through plots
for jplot=0,1 do begin

 ; define the margins
 xmar=[xmar1[jplot],xmar2[jplot]]

 ; define the and and y title
 xtitle='Evaporative fraction (observations)'
 if(jplot eq 0)then ytitle='Evaporative fraction (LSM simulations)' else ytitle=' '

 ; make a base plot
 plot, indgen(5), xrange=[0,xmax], yrange=[0,ymax], xstyle=1, ystyle=1, $
  xticks=3, xtitle=xtitle, ytitle=ytitle, xmargin=xmar, ymargin=ymar, /nodata
 plots, [0,xmax], [0,xmax]

 ; plot Budyko predictions
 oplot, eOverP_obs, budykoPred, psym=sym(1), symsize=2

 ; loop through realizations
 for iReal=0,nReal-2 do begin

  ; get the ratio of evaporation to precipitation
  eOverP_mod = reform((Qle[iReal  ,*]/LHvap)/Precip[*])

  ; define symbols and colors
  jxColor = iColor[ixType[iReal]]    ; color
  jxSymTy = iSymType[ixType[iReal]]  ; symbol type
  jxSymSz = iSymSize[ixType[iReal]]  ; symbol size

  ; plot model simulations
  if(jplot eq 0)then begin

   if(ixType[iReal] eq ixModel)then  begin
    ; plot the data
    oplot, eOverP_obs, eOverP_mod, psym=sym(model_symbols[iReal]), color=jxColor, symsize=2
    ; save the legend
    if(iReal eq 0)then legendCol=intarr(nModels)
    legendCol[iReal] = jxColor
    legendText = model_legend
    legendSym  = model_symbols
   endif

  ; plot physical and statistical benchmarks
  endif else begin
   if(ixType[iReal] ne ixModel)then begin
    ; plot the data
    oplot, eOverP_obs, eOverP_mod, psym=sym(benchmark_symbols[iReal-nModels]), color=jxColor, symsize=jxSymSz
    ; save the legend
    if(iReal eq nModels)then legendCol=intarr(nPhys+nStat)
    legendText = benchmark_legend
    legendSym  = benchmark_symbols
    legendCol[iReal-nModels] = jxColor 
   endif

  endelse

  ; save model simulations and obs
  if(jplot eq 0)then begin
   eOverP_save[iReal,  *] = eOverP_mod[*]
   eOverP_save[nReal-1,*] = eOverP_obs[*]
  endif

 endfor  ; looping through realizations

 ; plot a legend
 ytop = 2.9
 yInc = 0.1
 for iLegend=0,n_elements(legendText) do begin
  ypos = ytop - iLegend*yInc
  if(iLegend lt n_elements(legendText))then begin
   plots, 0.1, ypos, psym=sym(legendSym[iLegend]), color=legendCol[iLegend], symsize=2
   xyouts, 0.15, ypos-0.02, legendText[iLegend]
  endif else begin
   plots, 0.1, ypos, psym=sym(1), symsize=2
   xyouts, 0.15, ypos-0.02, 'Budyko curve'
  endelse
 endfor

endfor ; looping through plot types

; *****
; * MAKE BUDYKO BENCHMARK...
; **************************
window, 1, xs=1000, ys=1000, retain=2
device, decomposed=0
!P.MULTI=1

; define the rmse
rmsePenman = dblarr(nReal)

; compute the RMSE using all sites
rmsePenman[nReal-1] = sqrt(mean((budykoPred[*] - eOverP_save[nReal-1,*])^2.d))
for iReal=0,nReal-2 do begin
 rmsePenman[iReal] = sqrt(mean((eOverP_save[iReal,*] - eOverP_save[nReal-1,*])^2.d))
endfor

; re-sample
nTrial=1000
rmseResample = dblarr(nReal,nTrial)
for iTrial=0,nTrial-1 do begin
 ; get a random selection of sites
 jxSites = floor(randomu(seed,nSites) * float(nSites))
 ; compute the rmse using the random selection of sites
 for iReal=0,nReal-2 do begin
  rmseResample[iReal,iTrial] = sqrt(mean((eOverP_save[iReal,jxSites] - eOverP_save[nReal-1,jxSites])^2.d))
 endfor
 rmseResample[nReal-1,iTrial] = sqrt(mean((budykoPred[jxSites] - eOverP_save[nReal-1,jxSites])^2.d))
 ; get the resampled cdf
 ;yy = (dindgen(nReal-nPhys)+1.d)/double(nReal-nPhys)
 ;xx = [rmseResample[0:nModels-1],rmseResample[nModels+nPhys:nModels+nPhys+nStat-1],rmseResample[nReal-1]]
 ;ix = sort(xx)
 ;cdfResample[ix,iTrial] = yy
endfor  ; looping through trials

; get a delimiter
delim = replicate('-----', nReal)

; print header
print, '----------', delim, format='(a10,1x,25(a5,1x))'
print, 'Sim Type =', ixType, format='(a10,1x,25(i5,1x))' 
print, '----------', delim, format='(a10,1x,25(a5,1x))'

; print out a table
for iSite=0,nSites-1 do begin
 print, site_names[iSite], eOverP_save[*,iSite], budykoPred[iSite], format='(a10,1x,25(f5.3,1x))' 
endfor
print, '----------', delim, format='(a10,1x,25(a5,1x))'
print, 'rmsePenman', rmsePenman, format='(a10,1x,25(f5.3,1x))'

; make a base plot
plot, rmsePenman, xrange=[0,0.6], yrange=[0,1.05], xtitle='RMSE Evaporative fraction', ytitle='Cumulative probability', $
 xstyle=1, ystyle=1, /nodata

; define the legend text
legendText = [model_legend,benchmark_legend[nPhys:nPhys+nStat-1],'Budyko']
jxType = [ixType[0:nModels-1],ixType[nModels+nPhys:nModels+nPhys+nStat-1],ixType[nReal-1]]

; compute the cumulative probability
xx = [rmsePenman[0:nModels-1],rmsePenman[nModels+nPhys:nModels+nPhys+nStat-1],rmsePenman[nReal-1]]
ix = sort(xx)
xx = xx[ix]
yy = (dindgen(nReal-nPhys)+1.d)/double(nReal-nPhys)

; plot the cumulative probability
oplot, xx, yy, color = 80

jReal=0
; plot additinal information
ixColor=[250,0,80,0]
for iReal=0,nReal-1 do begin
 if(ixType[iReal] eq ixModel or ixType[iReal] eq ixStat or ixType[iReal] eq ixObs)then begin
  ; plot symbol
  jxColor=ixColor[jxType[ix[jReal]]]
  plots, xx[jReal], yy[jReal], color=jxColor, psym=sym(1), symsize=2
  ; plot uncertainty
  cdf = reform(rmseResample[ix[jReal],*])
  cdf = cdf[sort(cdf)]
  p05 = cdf[floor(0.05d*double(nTrial))]
  p95 = cdf[floor(0.95d*double(nTrial))]
  plots, [p05,p95], [yy[jReal], yy[jReal]]
  ; plot legend
  xyouts, xx[jReal]-0.01, yy[jReal]-0.005, legendText[ix[jReal]], alignment=1, color=jxColor
  print, legendText[ix[jReal]], xx[jReal], yy[jReal], jxType[ix[jReal]]
  jReal = jReal+1
 endif
endfor

stop
end
