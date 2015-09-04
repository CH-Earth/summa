pro plot_modelCompare

; define plotting parameters
window, 0, xs=1400, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=[0,3,2,0,0]

; define the date format
dummy = label_date(date_format=['%D-%M'])

; define the HRU
iHRU=0

; define the latent heat of vaporization
LHvap=2501000.d   ; J kg-1

; define the number of seconds in a day
secprday=86400.d

; define file paths
file_path = '/d1/mclark/PLUMBER_data/model_output/'

; define named variables to switch between canopy evaporation and canopy storage
ixStor=1            ; canopy storage
ixEvap=2            ; canopy evaporation
ixSWnet=3           ; net shortwave radiation
ixLWnet=4           ; net shortwave radiation
ixLatHeat=5         ; latent heat flux
ixSenHeat=6         ; latent heat flux
ixVegTemp=7         ; vegetation temperature
ixSoilMoist=8       ; root zone soil moisture
ixSoilStress=9      ; soil stress factor
ixStomatalResist=10 ; stomatal resistance 

; define desired variable
;ixVar=ixLatHeat
ixVar=ixStomatalResist

; define variables
case ixVar of

 ; =========================================================================================

 ; *** canopy storage
 ixStor: begin

  ; define variable range
  ymin=0
  ymax=0.5

  ; define canopy interception
  cVarNames=['CanopInt','CanopInt','CanopInt','CanopInt','SurfStor','scalarCanopyLiq']

  ; multiplier
  xMult=[1.d,1.d,1.d,1.d,1.d,1.d]

 end  ; canopy storage

 ; =========================================================================================

 ; *** canopy evaporation
 ixEvap: begin

  ; define variable range
  ymin=0
  ymax=250

  ; define canopy evaporation
  ;          kg/m2/s   kg/m2/s   mm/day  kg/m2/s  kg/m2/s     W/m2
  cVarNames=['ECanop','ECanop','ECanop','ECanop','ECanop','scalarLatHeatCanopyEvap']

  ; define multiplier
  xMult=[LHvap,LHvap,-LHvap/secprday,LHvap,LHvap,-1.d]

 end  ; canopy evaporation

 ; =========================================================================================

 ; *** net shortwave radiation
 ixSWnet: begin

  ; define variable range
  ymin=0
  ymax=1000

  ; define net sw radiation
  cVarNames=['SWnet','SWnet','SWnet','SWnet','SWnet','SWnet']

  ; define multiplier
  xMult=[1.d,1.d,1.d,1.d,1.d,1.d]

 end  ; net sw radiation

 ; =========================================================================================

 ; *** net longwave radiation
 ixLWnet: begin

  ; define variable range
  ymin=-200
  ymax=0

  ; define net sw radiation
  cVarNames=['LWnet','LWnet','LWnet','LWnet','LWnet','LWnet']

  ; define multiplier
  xMult=[1.d,1.d,1.d,1.d,1.d,1.d]

 end  ; net lw radiation

 ; =========================================================================================

 ; latent heat flux
 ixLatHeat: begin

  ; define variable range
  ymin=0
  ymax=500

  ; define latent heat flux
  cVarNames=['Qle','Qle','Qle','Qle','Qle','Qle']

  ; define multiplier
  xMult=[1.d,1.d,-1.d,1.d,1.d,-1.d]

 end

 ; =========================================================================================

 ; sensible heat flux
 ixSenHeat: begin

  ; define variable range
  ymin=0
  ymax=500

  ; define latent heat flux
  cVarNames=['Qh','Qh','Qh','Qh','Qh','Qh']

  ; define multiplier
  xMult=[1.d,1.d,-1.d,1.d,1.d,-1.d]

 end

 ; =========================================================================================

 ; vegetation temperature
 ixVegTemp: begin

  ; define variable range
  ymin=260
  ymax=320

  ; define latent heat flux
  cVarNames=['VegT','VegT','VegT','VegT','RadT','scalarCanopyTemp']

  ; define multiplier
  xMult=[1.d,1.d,1.d,1.d,1.d,1.d]

 end

 ; =========================================================================================

 ; root zone soil moisture
 ixSoilMoist: begin

  ; define variable range
  ymin=0
  ymax=10000000

  ; define root zone soil moisture
  cVarNames=['Qle','Qle','RootMoist','RootMoist','SoilMoist_tot','Qle']

  ; define multiplier
  xMult=[-1.d,-1.d,1.d,1.d,1.d,1.d]

 end

 ; =========================================================================================

 ; soil stress
 ixSoilStress: begin

  ; define variable range
  ymin=0
  ymax=1

  ; define root zone soil moisture
  cVarNames=['Qle','Qle','Qle','Qle','fsmc_pft','scalarTranspireLim']

  ; define multiplier
  xMult=[1.d,1.d,-1.d,1.d,1.d,1.d]

 end

 ; =========================================================================================

 ; stomatal resistance
 ixStomatalResist: begin

  ; define variable range
  ymin=0.d
  ymax=0.005

  ; define latent heat flux
  cVarNames=['Qle','Qle','Qle','Qle','Qle','scalarStomResistSunlit']
  ;cVarNames=['Qle','Qle','Qle','Qle','Qle','scalarStomResistShaded']
  ;cVarNames=['Qle','Qle','Qle','Qle','Qle','scalarVPair']

  ; define multiplier
  xMult=[1.d,1.d,-1.d,1.d,1.d,1.d]

 end

 ; =========================================================================================

 else: stop, 'cannot find desired variable'

endcase


; define the model names
model_names = ['CABLE.2.0',                  $
               'CABLE_2.0_SLI.vxh599_r553',  $
               'CHTESSEL',                   $
               'COLASSiB.2.0',               $
               ;'ISBA_SURFEX_3l.SURFEX7.3',   $
               ;'ISBA_SURFEX_dif.SURFEX7.3',  $
               ;'JULES.3.1',                  $
               'JULES3.1_altP',              $
               ;'Mosaic.1',                   $
               ;'NOAH.2.7.1',                 $
               ;'Noah.3.2',                   $
               ;'NOAH.3.3',                   $
               ;'ORCHIDEE.trunk_r1401',       $
               'SUMMA.1.0.exp.01.test']


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

; loop through sites
for iSite=0,nSites-1 do begin

 ; loop through the years
 for iYear=1994,2006 do begin

  ; define if we have done a plot
  done_plot=0

  ; loop through models
  for iModel=nModels-1,0,-1 do begin

   ; define the variable name
   cVarName = cVarNames[iModel]

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

    if(cVarName eq 'fsmc_pft')then begin
     xMax = max(mean(xvar, dimension=2), iPFT)
     xVar = reform(xVar[iPFT,*])
    endif

    if(cVarName eq 'scalarStomResistSunlit' or cVarName eq 'scalarStomResistShaded' or cVarName eq 'scalarVPair')then begin
     xVar = 1./xVar
     ;stop
    endif

   ; close the netcdf file
   ncdf_close, ncFileID

   ; check if within the date range
   caldat, djulian_mod[1], jm, jd, jyear_start
   caldat, djulian_mod[ntime_mod-2], jm, jd, jyear_end
   if(iYear lt jyear_start or iYear gt jyear_end)then continue

   ; define the title
   plotTitle=model_names[iModel] + '!C' + site_names[iSite] + ' (' + strtrim(iYear,2) + ')'

   ; make a hovmuller plot
   make_hovmuller, iYear, djulian_mod, xVar, cVarName, ymin, ymax, plotTitle

   ; define if we have done a plot
   done_plot=1

  endfor  ; looping through models

  ; write figure
  if(done_plot eq 1)then begin
   write_png, 'figures/hovmuller_'+cVarNames[5]+'_'+site_names[iSite]+'_'+strtrim(iYear,2)+'.png', tvrd(true=1)
   stop
  endif

 endfor ; looping through years

endfor  ; looping through sites

stop
end


; *****
; make a hovmuller diagram for the day-time...
; ********************************************

pro make_hovmuller, iyear_desire, dtime, varPlot, cVarName, vmin, vmax, plotTitle

; get the dates
xSmall = 1.d-6
caldat, dTime-xSmall, im, id, iyyy, ih, imin, asec

; identify the first year
iSubset=where(iyyy eq iyear_desire, nSubset)
;iSubset=where(iyyy eq iyear_desire and (im ge 1 and im le 4), nSubset)
;iSubset=where(iyyy eq iyear_desire and (im eq 3 or im eq 3) and (id ge 17 and id le 27), nSubset)

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
 polyfill, [ahour0,ahour1,ahour1,ahour0], [djulian-0.5, djulian-0.5, djulian+0.5, djulian+0.5], color=icolor

 ; print progress
 ;print, djulian, ih[itime], imin[itime], ahour1, varPlot[itime], icolor, format='(f30.10,1x,2(i4,1x),2(f9.3,1x),i4)'
 ;if(itime gt 500)then stop

endfor

end
