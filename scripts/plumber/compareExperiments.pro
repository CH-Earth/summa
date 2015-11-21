pro compareExperiments

; define plotting parameters
window, 0, xs=1200, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
;!P.MULTI=[0,3,3,0,0]
!P.MULTI=1

; define the path to the plumber data
plumber_path = '/Volumes/d1/mclark/PLUMBER_data/'

; define the path to the model output
model_path = plumber_path + 'model_output/'

; define the path to the site data
site_path =  plumber_path + 'site_data/flux/'

; define the benchmark path
benchmark_path = plumber_path + 'benchmark_data/'

; define the model names
model_names = ['CABLE.2.0',                  $
               'CABLE_2.0_SLI.vxh599_r553',  $
               ;'CHTESSEL',                   $
               'COLASSiB.2.0',               $
               'ISBA_SURFEX_3l.SURFEX7.3',   $
               'ISBA_SURFEX_dif.SURFEX7.3',  $
               'JULES.3.1',                  $
               'JULES3.1_altP',              $
               'Mosaic.1',                   $
               'NOAH.2.7.1',                 $
               'Noah.3.2',                   $
               'NOAH.3.3',                   $
               ;'ORCHIDEE.trunk_r1401',       $
               'SUMMA.1.0.exp.02.001',       $
               'SUMMA.1.0.exp.02.002',       $
               'SUMMA.1.0.exp.02.003',       $
               'SUMMA.1.0.exp.02.004',       $
               'SUMMA.1.0.exp.02.005',       $
               'SUMMA.1.0.exp.02.006',       $
               'SUMMA.1.0.exp.02.007',       $
               'SUMMA.1.0.exp.02.008',       $
               'SUMMA.1.0.exp.02.009',       $
               'SUMMA.1.0.exp.02.010',       $
               'SUMMA.1.0.exp.02.011',       $
               'SUMMA.1.0.exp.02.012',       $
               'SUMMA.1.0.exp.02.013',       $
               'SUMMA.1.0.exp.02.014',       $
               'SUMMA.1.0.exp.02.015',       $
               'SUMMA.1.0.exp.02.016',       $
               'SUMMA.1.0.exp.02.017',       $
               'SUMMA.1.0.exp.02.018',       $
               'SUMMA.1.0.exp.02.019',       $
               'SUMMA.1.0.exp.02.021',       $
               'SUMMA.1.0.exp.02.022',       $
               'SUMMA.1.0.exp.02.030',       $
               'SUMMA.1.0.exp.02.031',       $
               'SUMMA.1.0.exp.02.032',       $
               'SUMMA.1.0.exp.02.000'       ]

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
nModels = n_elements(model_names)
nSites  = n_elements(site_names)

; define the number of stats
nStats=9

; define statistics
statName = ['Absolute bias', $
            'Variance ratio', $
            '1 - Correlation',    $
            'Normalized mean absolute error', $
            'Differnce in 5th percentile', $
            'Differnce in 95th percentile', $
            'Absolute difference in Skewness', $
            'Absolute difference in Kurtosis', $
            '1 - (Overlap statistic)']

; define labels for the statistics (for the filename)
statLabel = ['AbsBias', $
             'variRat', $
             '1minusR',    $
             'NormMAE', $
             'perc_05', $
             'perc_95', $
             'skewRat', $
             'kurtRat', $
             'overlap']

; define model labels
model_labels = strarr(nModels)
for iModel=0,nModels-1 do model_labels[iModel] = model_names[iModel] + '(' + strtrim(iModel,2) + ')'

; loop through sites
for iSite=0,nSites-1 do begin

 ; skip computations (done already)
 goto, got_stats

 ; *****
 ; * READ IN THE SITE OBS...
 ; *************************

 ; define filename
 site_filename = site_path + site_names[iSite] + 'Fluxnet.1.4_flux.nc'

 ; open file for reading
 nc_file = ncdf_open(site_filename, /nowrite)

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
  if(strtrim(tunit_words[0],2) eq 'seconds') then aoff=1.d/86400.d else stop, 'unknown time units'

  ; extract the time vector
  ncdf_varget, nc_file, ivar_id, atime
  djulian = bjulian + atime*aoff

  ; get the number of time elements
  ntime = n_elements(djulian)

  ; read in sensible heat
  ivar_id = ncdf_varid(nc_file,'Qh')
  ncdf_varget, nc_file, ivar_id, Qh
  Qh_obs = reform(Qh)

  ; read in latent heat
  ivar_id = ncdf_varid(nc_file,'Qle')
  ncdf_varget, nc_file, ivar_id, Qle
  Qle_obs = reform(Qle)

 ; close netcdf file
 ncdf_close, nc_file

 ; *****
 ; * READ IN THE MODEL SIMULATIONS...
 ; **********************************

 ; define array holding model simulations
 Qh_model = fltarr(nModels,ntime)
 Qle_model = fltarr(nModels,ntime)

 ; loop through models
 for imodel=0,nmodels-1 do begin

  ; define the model filename
  model_filename = model_path + model_names[imodel] + '/' + model_names[imodel] + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'
  print, imodel, ' - ', model_filename

  ; open file for reading
  nc_file = ncdf_open(model_filename, /nowrite)

   ; read in the sensible heat flux
   ivar_id = ncdf_varid(nc_file,'Qh')
   ncdf_varget, nc_file, ivar_id, Qh
   Qh_model[imodel,*] = reform(Qh)

   ; read in the latent heat flux
   ivar_id = ncdf_varid(nc_file,'Qle')
   ncdf_varget, nc_file, ivar_id, Qle
   Qle_model[imodel,*] = reform(Qle)

  ; close netcdf file
  ncdf_close, nc_file 

  ; change sign for CHTESSEL and SUMMA
  if(model_names[imodel] eq 'CHTESSEL' or strmid(model_names[imodel],0,5) eq 'SUMMA')then begin
   Qh_model[imodel,*]  = -Qh_model[imodel,*]
   Qle_model[imodel,*] = -Qle_model[imodel,*]
  endif

 endfor  ; (looping through models)

 ; *****
 ; * COMPUTE THE STATISTICS...
 ; ***************************

 ; define an array for the statistics
 xStats_Qh = replicate(0.d, nModels,nModels,nstats)
 xStats_Qle = replicate(0.d, nModels,nModels,nstats)

 ; loop through models
 for iModel=0,nModels-1 do begin
  for jModel=0,iModel do begin

   ; compute stats (compare models against obs)
   if(imodel eq jmodel)then begin
    compute_stats, reform(Qh_model[iModel,*]),  Qh_obs,  Qh_statsModel
    compute_stats, reform(Qle_model[iModel,*]), Qle_obs, Qle_statsModel

   ; compute statistics (compare models against models)
   endif else begin
    compute_stats, reform(Qh_model[iModel,*]),  reform(Qh_model[jModel,*]),  Qh_statsModel    
    compute_stats, reform(Qle_model[iModel,*]), reform(Qle_model[jModel,*]), Qle_statsModel  
   endelse

   ; save statistics
   xStats_Qh[iModel,jModel,*] = Qh_statsModel[*]
   xStats_Qle[iModel,jModel,*] = Qle_statsModel[*]

  endfor  ; looping through models

  ; print progress
  print, 'computing stats for model ', iModel, max(xStats_Qle[iModel,*,0])

 endfor  ; looping through models

 ; save statistics
 save, /variables, filename='xIDLsave/modelStatistics_'+site_names[iSite]+'.sav'

 ; label if have the statistics already
 got_stats:

 ; restore the statistics
 restore, 'xIDLsave/modelStatistics_'+site_names[iSite]+'.sav'

 ; *****
 ; * PLOT THE STATISTICS...
 ; ************************

 ; define plot range for each statistic
 vmin = [ 0, 0, 0, 0,  0,   0, 0, 0, 0]
 vmax = [50, 2, 1, 2, 75, 300, 5, 5, 1]

 ; loop through the statistics
 for iStat=0,nStats-1 do begin

  ; define the range
  xMin = vMin[iStat]
  xMax = vMax[iStat]

  ; make a base plot
  plot, indgen(nModels), xrange=[0,nModels], yrange=[0,nModels], xstyle=1, ystyle=1, $
   ytickname=[model_labels,' '], yticks=nModels, xmargin=[30,2], xticklen=(-0.02), $
   title=site_names[iSite] + ' [' + statName[iStat] + ']', /nodata

  ; plot the results
  for iModel=0,nModels-1 do begin
   for jModel=0,iModel do begin
    ; define the color
    icolor= 50. + 250. * (xStats_Qle[iModel,jModel,iStat] - xMin) / (xMax - xMin)
    if(xStats_Qle[iModel,jModel,iStat] lt xMin)then iColor=255
    if(xStats_Qle[iModel,jModel,iStat] gt xMax)then iColor=255
    ; plot the data
    xx = [iModel,iModel+1,iModel+1,iModel]
    yy = [jModel,jModel,jModel+1,jModel+1]
    polyfill, xx, yy, color=icolor
   endfor
  endfor

  ; define name of figure
  figname = 'plumber_modelComparison_' + statLabel[iStat] + '_' + site_names[iSite] + '.png'

  ; write figure
  write_png, 'figures/modelComparison/'+figname, tvrd(true=1)

  ; test
  ;if(statLabel[istat] eq 'kurtRat')then stop

 endfor

endfor  ; looping through sites


stop
end


; *****
; * COMPUTE THE STATISTICS...
; ***************************
pro compute_stats, sim, obs, stats

; define the statistics vector
stats = dblarr(9)

; compute the first few moments
xx = moment(sim, mean=simMean, sdev=simSdev, skewness=simSkew, kurtosis=simKurt, mdev=simMdev, /double)
xx = moment(obs, mean=obsMean, sdev=obsSdev, skewness=obsSkew, kurtosis=obsKurt, mdev=obsMdev, /double)

; sort the data
simSort = sim[sort(sim)]
obsSort = obs[sort(obs)]

; get the indices for the 5 and 95 percentiles
num = n_elements(sim)
i05 = floor(double(num)*0.05d)
i95 = floor(double(num)*0.95d)

; compute the overlap statistic
nBin = 25
xOvr = 0.d ; initialize the overlap statistic
xMin = min([min(obs),min(sim)])
xMax = max([max(obs),max(sim)])
xInc = (xMax - xMin)/double(nBin)
for iBin=0,nBin-1 do begin
 xxBin = xMin + double(iBin)*xInc
 ixBin = where(sim ge xxBin and sim lt xxBin+xInc, nSim)
 ixBin = where(obs ge xxBin and obs lt xxBin+xInc, nObs)
 xFreq = min([double(nSim)/double(num),double(nObs)/double(num)])
 xOvr  = xOvr + xfreq
endfor

; compute the statistics
stats[0] = abs(simMean - obsMean)            ; absolute bias
stats[1] = abs(1.d - simSdev/obsSdev)        ; variance ratio
stats[2] = 1.d - correlate(sim,obs)          ; correlation
stats[3] = mean(abs(obs-sim))/obsMdev        ; normalized mean error
stats[4] = abs(simSort[i05] - obsSort[i05])  ; 5% error
stats[5] = abs(simSort[i95] - obsSort[i95])  ; 95% error
stats[6] = abs(1.d - simSkew/obsSkew)        ; skewness
stats[7] = abs(1.d - simKurt/obsKurt)        ; kurtosis
stats[8] = 1.d - xOvr                        ; overlap

end
