pro plot_skillStatistics

; define plotting parameters
window, 0, xs=2000, ys=1400, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3
!P.COLOR=0
erase, color=255
!P.MULTI=[0,3,3,0,0]

; define the path to the plumber data
plumber_path = '/d1/mclark/PLUMBER_data/'

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
               'SUMMA.1.0.exp.02.000',       $
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
               'SUMMA.1.0.exp.02.012'       ]

; define the benchmark names
sb_names = ['1lin','2lin','3km27']  ; statistical benchmarks
pb_names = ['Manabe_Bucket.2','Penman_Monteith.1'] ; physical benchmarks

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

; define the number of statistical and physical benchmarks
nSB     = n_elements(sb_names)  ; statistical benchmarks
nPB     = n_elements(pb_names)  ; physical benchmarks
nBench  = nSB + nPB

; define arrays for the stats
nStats=9
Qh_stats = dblarr(nSites,nBench+nModels,nStats)
Qle_stats = dblarr(nSites,nBench+nModels,nStats)

; skip computations (done already)
;goto, got_stats

; loop through sites
for iSite=0,nSites-1 do begin

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

  ; read in sensible heat quality control flag
  ivar_id = ncdf_varid(nc_file,'Qh_qc')
  ncdf_varget, nc_file, ivar_id, Qh_qc
  Qh_qc = floor(reform(Qh_qc)+0.5)

  ; read in latent heat
  ivar_id = ncdf_varid(nc_file,'Qle')
  ncdf_varget, nc_file, ivar_id, Qle
  Qle_obs = reform(Qle)

  ; read in latent heat quality control flag
  ivar_id = ncdf_varid(nc_file,'Qle_qc')
  ncdf_varget, nc_file, ivar_id, Qle_qc
  Qle_qc = floor(reform(Qle_qc)+0.5)

 ; close netcdf file
 ncdf_close, nc_file

 ; *****
 ; * READ IN THE BENCHMARKS...
 ; ***************************

 ; initialize benchmark index
 iBench=0

 ; define array holding benchmarks
 Qh_bench = fltarr(nBench,ntime)
 Qle_bench = fltarr(nBench,ntime)

 ; define filename
 benchmark_filename = benchmark_path + site_names[iSite] + 'Fluxnet_1.4_PLUMBER_benchmarks.nc'

 ; open file for reading
 nc_file = ncdf_open(benchmark_filename, /nowrite)

  ; loop thru statistical benchmarks
  for iStat=0,nSB-1 do begin

   ; read in sensible heat
   ivar_id = ncdf_varid(nc_file,'Qh_'+sb_names[iStat])
   ncdf_varget, nc_file, ivar_id, Qh
   Qh_bench[iBench,*] = reform(Qh)

   ; read in latent heat
   ivar_id = ncdf_varid(nc_file,'Qle_'+sb_names[iStat])
   ncdf_varget, nc_file, ivar_id, Qle
   Qle_bench[iBench,*] = reform(Qle)

   ; increment benchmark index
   iBench = iBench + 1

  endfor  ; looping through statistical benchmarks

 ; close netcdf file
 ncdf_close, nc_file

 ; loop through the physical benchmarks
 for imodel=0,nPB-1 do begin

  ; define the model filename
  benchmark_filename = model_path + pb_names[imodel] + '/' + pb_names[imodel] + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'
  print, benchmark_filename

  ; open file for reading
  nc_file = ncdf_open(benchmark_filename, /nowrite)

   ; read in the sensible heat flux
   ivar_id = ncdf_varid(nc_file,'Qh')
   ncdf_varget, nc_file, ivar_id, Qh
   Qh_bench[iBench,*] = reform(Qh)

   ; read in the latent heat flux
   ivar_id = ncdf_varid(nc_file,'Qle')
   ncdf_varget, nc_file, ivar_id, Qle
   Qle_bench[iBench,*] = reform(Qle)

  ; close netcdf file
  ncdf_close, nc_file

  ; increment benchmark index
  iBench = iBench + 1

 endfor  ; looping through the physical benchmarks

 ; *****
 ; * GET THE MODEL SIMULATIONS...
 ; ******************************

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

 ; compute statistics for the benchmarks
 for iBench=0,nBench-1 do begin
  compute_stats, reform(Qh_bench[iBench,*]), Qh_obs, Qh_statsBench     & Qh_stats[iSite,iBench,*] = Qh_statsBench[*]
  compute_stats, reform(Qle_bench[iBench,*]), Qle_obs, Qle_statsBench  & Qle_stats[iSite,iBench,*] = Qle_statsBench[*]
 endfor  ; looping through benchmarks

 ; compute the statistics for the models
 for iModel=nModels-1,0,-1 do begin
  compute_stats, reform(Qh_model[iModel,*]), Qh_obs, Qh_statsModel     & Qh_stats[iSite,nBench+iModel,*] = Qh_statsModel[*]
  compute_stats, reform(Qle_model[iModel,*]), Qle_obs, Qle_statsModel  & Qle_stats[iSite,nBench+iModel,*] = Qle_statsModel[*]
 endfor

endfor  ; looping through sites

save, /variables, filename='xIDLsave/skillStatistics.sav'

got_stats:
restore, 'xIDLsave/skillStatistics.sav'

; *****
; * PLOT THE STATISTICS...
; ************************

; skip the plot
;goto, got_plot

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

; define plot range for each statistic
vmin = [ 0, 0, 0, 0,  0,   0,  0,   0, 0]
vmax = [50, 2, 1, 2, 75, 300, 10, 150, 1]

; loop through the statistics
for iStat=0,nStats-1 do begin

 ; make a base plot (for a given site)
 plot, indgen(nSites+2), xrange=[0,nSites+1], yrange=[vmin[iStat],vmax[iStat]], xstyle=9, ystyle=1, xticklen=(-0.02),$
  xtickname=[' ',strtrim(indgen(nSites)+1,2),' '], xticks=nSites+1, ytitle=statName[iStat], $
  /nodata
 plots, [0,nSites+1], [vmax[iStat],vmax[iStat]]

 ; loop through sites
 for iSite=0,nSites-1 do begin

  ; plot the model simulations
  for iModel=0,nModels-1 do begin
   plots, float(iSite+1)-0.25,  Qh_stats[iSite,nBench+iModel,iStat], psym=sym(1), color=250
   plots, float(iSite+1)+0.25,  Qle_stats[iSite,nBench+iModel,iStat], psym=sym(1), color=80
  endfor

  ; plot the benchmarks
  for iBench=0,nBench-1 do begin
   plots, float(iSite+1)-0.25,  Qh_stats[iSite,iBench,iStat], psym=sym(iBench+6), color=250, symsize=2
   plots, float(iSite+1)+0.25,  Qle_stats[iSite,iBench,iStat], psym=sym(iBench+6), color=80, symsize=2
  endfor

  ; print out results for Loobos
  if(site_names[iSite] eq 'Loobos')then begin
   print, 'SUMMA: ', statName[iStat], Qh_stats[iSite,nBench+nModels-1,iStat], Qle_stats[iSite,nBench+nModels-1,iStat], $
    format='(a10,1x,a35,1x,2(f15.10,1x))'
   print, '1lin: ', statName[iStat], Qh_stats[iSite,0,iStat], Qle_stats[iSite,0,iStat], $
    format='(a10,1x,a35,1x,2(f15.10,1x))'
   print, '2lin: ', statName[iStat], Qh_stats[iSite,1,iStat], Qle_stats[iSite,1,iStat], $
    format='(a10,1x,a35,1x,2(f15.10,1x))'
   print, '3km27: ', statName[iStat], Qh_stats[iSite,2,iStat], Qle_stats[iSite,2,iStat], $
    format='(a10,1x,a35,1x,2(f15.10,1x))'
   print, 'MANABE: ', statName[iStat], Qh_stats[iSite,3,iStat], Qle_stats[iSite,3,iStat], $
    format='(a10,1x,a35,1x,2(f15.10,1x))'
  endif

 endfor  ; looping through the sites

endfor  ; looping through the statistics

; write figure
write_png, 'figures/plumber_skill.png', tvrd(true=1)

; jump point
got_plot:

; *****
; * PLOT THE RANKING...
; *********************

; define the ranking
Qh_rank = dblarr(nBench+1,nModels)
Qle_rank = dblarr(nBench+1,nModels)

; define the desired statistics
ixStat=[0,1,2,3]
nxStat=n_elements(ixStat)

; loop through the models
for iModel=0,nModels-1 do begin

 ; initialize Qh_rank and Qle_rank
 Qh_rank[*,iModel] = 0.d
 Qle_rank[*,iModel] = 0.d

 ; loop through the sites and statistics
 for iSite=0,nSites-1 do begin
  for jStat=0,nxStat-1 do begin

    ; define the index of the statistic
    iStat=ixStat[jstat]

    ; get a vector of six values
    Qh_vector =  [reform(Qh_stats[iSite,0:nBench-1,iStat]),Qh_stats[iSite,nBench+iModel,iStat]]
    Qle_vector = [reform(Qle_stats[iSite,0:nBench-1,iStat]),Qle_stats[iSite,nBench+iModel,iStat]]

    ; sort the six values
    Qh_sort = sort(Qh_vector)
    Qle_sort = sort(Qle_vector)

    ; get the ranking
    Qh_vector[Qh_sort]   = dindgen(nBench+1)+1.d
    Qle_vector[Qle_sort] = dindgen(nBench+1)+1.d

    ; accumulate the ranking
    Qh_rank[*,iModel] = Qh_rank[*,iModel] + Qh_vector[*]
    Qle_rank[*,iModel] = Qle_rank[*,iModel] + Qle_vector[*]

    ; print summa
    ;if(model_names[iModel] eq 'SUMMA.1.0')then begin
    ; if(jstat eq 0)then print, '**'
    ; print, 'Qh rank: ', site_names[isite], istat, Qh_vector[*] 
    ; print, 'Qle rank: ', site_names[isite], istat, Qle_vector[*]
    ;endif

  endfor  ; (looping through stats)
 endfor  ; (looping through sites)

 ; normalize
 Qh_rank[*,iModel] = Qh_rank[*,iModel] / double(nSites*nxStat)
 Qle_rank[*,iModel] = Qle_rank[*,iModel] / double(nSites*nxStat)

 ; print results
 print, '**'
 print, model_names[iModel], Qh_rank[*,iModel], format='(a25,1x,6(f9.3,1x))'
 print, model_names[iModel], Qle_rank[*,iModel], format='(a25,1x,6(f9.3,1x))'


endfor  ; (loop through models)


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
