pro plot_hovmuller

; define plotting parameters
window, 0, xs=2000, ys=1400, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=[0,5,4,0,0]

; define the date format
dummy = label_date(date_format=['%D-%M'])

; define the HRU
iHRU=0

; define variable
;cVarName='Qle'
cVarName='scalarLatHeatCanopyEvap'
;cVarName='scalarLatHeatCanopyTrans'

; define variable range
ymin=-500
ymax=500

; define the model names
model_names = ['CABLE.2.0',                  $
               ;'CABLE_2.0_SLI.vxh599_r553',  $
               ;'CHTESSEL',                   $
               ;'COLASSiB.2.0',               $
               ;'ISBA_SURFEX_3l.SURFEX7.3',   $
               ;'ISBA_SURFEX_dif.SURFEX7.3',  $
               ;'JULES.3.1',                  $
               ;'JULES3.1_altP',              $
               ;'Mosaic.1',                   $
               ;'NOAH.2.7.1',                 $
               ;'Noah.3.2',                   $
               ;'NOAH.3.3',                   $
               ;'ORCHIDEE.trunk_r1401',       $
               'SUMMA.1.0.exp.01.003']

; flux multiplier
xMult=[1.d,-1.d]

; define the site names
site_names = ['Amplero',     $
              'Blodgett',    $
              'Bugac',       $
              'ElSaler',     $
              'ElSaler2',    $
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

site_pfts = ['Grassland',          $        ; 'Amplero',     
             'Evergreen Needleleaf',$       ; 'Blodgett',    
             'Grassland',          $        ; 'Bugac',       
             'Evergreen Needleleaf',$       ; 'ElSaler',    
             'Cropland',           $        ; 'ElSaler2',     
             'Evergreen Broadleaf', $       ; 'Espirra',     
             'Grassland',          $        ; 'FortPeck',    
             'Deciduous Broadleaf', $       ; 'Harvard',     
             'Deciduous Broadleaf', $       ; 'Hesse',       
             'Woody Savanna',       $       ; 'Howard',      
             'Evergreen Needleleaf',$       ; 'Howlandm',    
             'Evergreen Needleleaf',$       ; 'Hyytiala',    
             'Savanna',            $        ; 'Kruger',      
             'Evergreen Needleleaf',$       ; 'Loobos',      
             'Permanent Wetland',   $       ; 'Merbleue',    
             'Woody Savanna',       $       ; 'Mopane',      
             'Evergreen Broadleaf', $       ; 'Palang',      
             'Mixed Forest',        $       ; 'Sylvania',    
             'Evergreen Broadleaf', $       ; 'Tumba',       
             'Deciduous Broadleaf'  ]       ; 'UniMich'      

; define the number of models and sites
nModels = n_elements(model_names)
nSites  = n_elements(site_names)

; define file paths
file_path = '/d1/mclark/PLUMBER_data/model_output/'

; loop through models
for iModel=nModels-1,0,-1 do begin

 ; loop through sites
 for iSite=0,nSites-1 do begin

  ; define the file name
  file_name = model_names[iModel] + '/' + model_names[iModel] + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'

  ; open files
  ncFileID = ncdf_open(file_path+file_name, /nowrite)

   ; get time units
   ivar_id = ncdf_varid(ncFileID,'time')
   ncdf_attget, ncFileID, ivar_id, 'units', bunits
   cunits = string(bunits)

   ; extract the units "words"
   tunit_words = strsplit(string(cunits),' ',/extract)
   tunit_idate = fix(strsplit(tunit_words[2],'-',/extract))
   tunit_ihour = fix(strsplit(tunit_words[3],':',/extract))
   bjulian     = julday(tunit_idate[1],tunit_idate[2],tunit_idate[0],tunit_ihour[0],tunit_ihour[1],tunit_ihour[2])

   ; get the offset in days
   if(strtrim(tunit_words[0],2) eq 'seconds') then aoff=1.d/86400.d else stop, 'unknown time units'

   ; extract the time vector
   ncdf_varget, ncFileID, ivar_id, atime
   djulian_mod = bjulian + atime*aoff

   ; get the number of time elements
   ntime_mod = n_elements(djulian_mod)

   ; get the desired variable
   ivar_id = ncdf_varid(ncFileID,cVarName)
   ncdf_varget, ncFileID, ivar_id, xVar

   ; modify the variable
   xVar=reform(xVar)*xMult[iModel]

  ; close the netcdf file
  ncdf_close, ncFileID

  ; define the year
  caldat, djulian_mod[1], jm, jd, jyyy

  ; define the title
  plotTitle=model_names[iModel] + '!C' + site_names[iSite] + ' (' + strtrim(jyyy,2) + ')'

  ; make a hovmuller plot
  make_hovmuller, djulian_mod, xVar, cVarName, ymin, ymax, plotTitle

 endfor  ; looping through sites

 ; write figure
 write_png, 'figures/hovmuller_'+cVarName+'_'+model_names[iModel]+'.png', tvrd(true=1)

 stop
endfor  ; looping through models

stop
end


; *****
; make a hovmuller diagram for the day-time...
; ********************************************

pro make_hovmuller, dtime, varPlot, cVarName, vmin, vmax, plotTitle

; get the dates
xSmall = 1.d-6
caldat, dTime-xSmall, im, id, iyyy, ih, imin, asec

; identify the first year
iSubset=where(iyyy eq iyyy[0], nSubset)

; define the simulation time
time0 = dTime[iSubset[0]]
time1 = dTime[iSubset[nSubset-1]]

; define the x-tick labels
xTickLabels = [' ',[strtrim(indgen(7)*3+3,2)],' ']


; make a base plot
plot, indgen(24)+1, xrange=[0,24], yrange=[time0,time1], xstyle=9, ystyle=1, $
 xmargin = [16,2], xticks=8, xtickname=xTickLabels, ytitle=cVarName, $
 xcharsize=1.5, ycharsize=1.5, xticklen=(-0.02), title=plotTitle, $
 ytickformat=['label_date'], yticks=5, yticklen=(-0.02), $
 ymargin=[3,5], /nodata
plots, [0,24], [time1,time1]

; loop through the time series
for jTime=0,nSubset-1 do begin

 ; get the time index
 iTime = iSubset[jTime]

 ; get the julian day
 djulian = julday(im[itime], id[itime], iyyy[itime])

 ; identify the color
 if(varPlot[itime] gt -9998.)then begin
  icolor = ( (varPlot[itime] - vmin) / (vmax - vmin) )*200.d + 50.d
  if(varPlot[itime] lt vmin)then icolor=50
  if(varPlot[itime] gt vmax)then icolor=250
 endif else begin
  icolor = 255
 endelse

 ; get the decimal time
 aHour1 = double(ih[itime]) + double(imin[itime])/60.d + asec[itime]/3600.d
 aHour0 = aHour1 - 0.5d

 ; plot data
 plots, [ahour0,ahour1,ahour1,ahour0], [djulian-0.5, djulian-0.5, djulian+0.5, djulian+0.5], color=icolor

 ; print progress
 ;print, djulian, ih[itime], imin[itime], ahour1, varPlot[itime], icolor, format='(f30.10,1x,2(i4,1x),2(f9.3,1x),i4)'
 ;if(itime gt 500)then stop

endfor


end
