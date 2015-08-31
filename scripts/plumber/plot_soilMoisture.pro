pro plot_soilMoisture

; define plotting parameters
window, 0, xs=1400, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=[0,1,2,0,0]

; define the date format
dummy = label_date(date_format=['%D-%M!C%Y'])

; refine colors
tvlct, r, g, b, /get
r = reverse(r)
g = reverse(g)
b = reverse(b)
r[0] = 0
g[0] = 0
b[0] = 0
r[255] = 255
g[255] = 255
b[255] = 255
tvlct, r, g, b

; define the HRU
iHRU=0

; define the file path
file_path = '/d1/mclark/PLUMBER_data/model_output/'

; define the model names
cable_name = 'CABLE_2.0_SLI.vxh599_r553'
summa_name = 'SUMMA.1.0.exp.01.test'

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

site_pfts = ['Grassland',           $       ; 'Amplero',     
             'Evergreen Needleleaf',$       ; 'Blodgett',    
             'Grassland',           $       ; 'Bugac',       
             'Evergreen Needleleaf',$       ; 'ElSaler',    
             'Cropland',            $       ; 'ElSaler2',     
             'Evergreen Broadleaf', $       ; 'Espirra',     
             'Grassland',           $       ; 'FortPeck',    
             'Deciduous Broadleaf', $       ; 'Harvard',     
             'Deciduous Broadleaf', $       ; 'Hesse',       
             'Woody Savanna',       $       ; 'Howard',      
             'Evergreen Needleleaf',$       ; 'Howlandm',    
             'Evergreen Needleleaf',$       ; 'Hyytiala',    
             'Savanna',             $       ; 'Kruger',      
             'Evergreen Needleleaf',$       ; 'Loobos',      
             'Permanent Wetland',   $       ; 'Merbleue',    
             'Woody Savanna',       $       ; 'Mopane',      
             'Evergreen Broadleaf', $       ; 'Palang',      
             'Mixed Forest',        $       ; 'Sylvania',    
             'Evergreen Broadleaf', $       ; 'Tumba',       
             'Deciduous Broadleaf'  ]       ; 'UniMich'      

; define the number of models and sites
nSites  = n_elements(site_names)

; loop through sites
for iSite=0,nSites-1 do begin


 ; ***********
 ; ** CABLE...
 ; ***********

 ; define the file name
 file_name = cable_name + '/' + cable_name + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'

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
  djulian_cable = bjulian + atime*aoff

  ; get the soil moisture
  ivar_id = ncdf_varid(ncFileID,'SoilMoist')
  ncdf_varget, ncFileID, ivar_id, xVar
  soilMoist=reform(xVar)

  ; get the depth of each soil layer
  ivar_id = ncdf_varid(ncFileID,'zse')
  ncdf_varget, ncFileID, ivar_id, xVar
  soilDepth=reform(xVar)

 ; close the netcdf file
 ncdf_close, ncFileID

 ; define the soil height
 nSoil = n_elements(soilDepth)
 soilHeight = replicate(0.d, nSoil+1)
 for iSoil=0,nSoil-1 do soilHeight[iSoil+1] = soilHeight[iSoil] + soilDepth[iSoil]
 print, soilHeight
 
 ; define the title
 plotTitle='CABLE: ' + site_names[iSite] 

 ; identify the time period to plot
 iSubset=where(djulian_cable gt 0, nSubset)

 ; define the simulation time
 time0 = djulian_cable[iSubset[0]]
 time1 = djulian_cable[iSubset[nSubset-1]]

 ; make a base plot
 plot, indgen(5), xrange=[time0,time1], yrange=[5,0], xstyle=9, ystyle=1, $
  xtickformat=['label_date'], xticks=9, xticklen=(-0.02), $
  title=plotTitle, /nodata
 plots, [time1,time1], [5,0]

 ; define the plot range
 vmin = 0.0
 vmax = 0.5

 ; make a hovmuller plot
 time_hovmuller, djulian_cable, soilHeight, soilMoist, vmin, vmax, plotTitle

 ; ***********
 ; ** SUMMA...
 ; ***********

 ; define the title
 plotTitle='SUMMA: ' + site_names[iSite] 

 ; make a base plot
 plot, indgen(5), xrange=[time0,time1], yrange=[5,0], xstyle=9, ystyle=1, $
  xtickformat=['label_date'], xticks=9, xticklen=(-0.02), $
  title=plotTitle, /nodata
 plots, [time1,time1], [5,0]

 ; get the list of files
 spawn, 'ls -1 ' + file_path + summa_name + '/orig/*' + site_names[iSite] + '*initialPlumberTest.nc', file_list

 ; loop through files
 for iFile=0,n_elements(file_list)-1 do begin

  ; print progress
  print, file_list[ifile]

  ; open files
  ncFileID = ncdf_open(file_list[ifile], /nowrite)

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
   djulian_summa = bjulian + atime*aoff

   ; extract the number of snow layers
   ivar_id = ncdf_varid(ncFileID,'nSnow')
   ncdf_varget, ncFileID, ivar_id, xVar
   nSnow=reform(xVar)

   ; extract the total number of layers
   ivar_id = ncdf_varid(ncFileID,'nLayers')
   ncdf_varget, ncFileID, ivar_id, xVar
   nLayers=reform(xVar)

   ; extract the start index for the layer midpoint
   ivar_id = ncdf_varid(ncFileID,'midTotoStartIndex')
   ncdf_varget, ncFileID, ivar_id, xVar
   midTotoStartIndex=reform(xVar)

   ; extract the start index for the layer interfaces
   ivar_id = ncdf_varid(ncFileID,'ifcTotoStartIndex')
   ncdf_varget, ncFileID, ivar_id, xVar
   ifcTotoStartIndex=reform(xVar)

   ; get the number of soil layers and time
   nSoil = min(nLayers) - min(nSnow)
   nTime = n_elements(djulian_summa)

   ; get the matrix of soil moisture
   soilMoist = fltarr(nSoil,nTime)

   ; loop through time
   for iTime=0,nTime-1 do begin

    ; get the height
    if(iTime eq 0)then begin
     ivar_id = ncdf_varid(ncFileID,'iLayerHeight')
     ncdf_varget, ncFileID, ivar_id, xVar, offset=[iHRU,ifcTotoStartIndex[itime]+nSnow[iTime]-1], count=[1,nSoil+1]
     soilHeight = reform(xVar)
    endif

    ; get the soil moisture
    ivar_id = ncdf_varid(ncFileID,'mLayerVolFracLiq')
    ncdf_varget, ncFileID, ivar_id, xVar, offset=[iHRU,midTotoStartIndex[itime]+nSnow[iTime]-1], count=[1,nSoil]
    soilMoist[*,iTime] = xVar[0,*]

   endfor  ; looping through time

   ; make a hovmuller plot
   time_hovmuller, djulian_summa, soilHeight, soilMoist, vmin, vmax, plotTitle

 endfor   ; looping through files

 ; write figure
 write_png, 'figures/soilMoist_'+site_names[iSite]+'.png', tvrd(true=1)
 stop

endfor  ; looping through sites

stop
end


; *****
; make a hovmuller diagram for the day-time...
; ********************************************

pro time_hovmuller, dTime, soilHeight, soilMoist, vmin, vmax, plotTitle

; get the dates
xSmall = 1.d-6
caldat, dTime-xSmall, im, id, iyyy, ih, imin, asec

; loop through the time series
for iTime=0,n_elements(dTime)-1 do begin

 ; get the julian day
 djulian = julday(im[itime], id[itime], iyyy[itime], ih[itime], imin[itime])

 ; loop through the soil depths
 for iDepth=0,n_elements(soilHeight)-2 do begin  ; -2 because height is at the interface

  ; identify the color
  if(soilMoist[iDepth,iTime] gt -9998.)then begin
   icolor = ( (soilMoist[iDepth,iTime] - vmin) / (vmax - vmin) )*200.d + 50.d
   if(soilMoist[iDepth,iTime] lt vmin)then icolor=50
   if(soilMoist[iDepth,iTime] gt vmax)then icolor=250
   ;if(soilMoist[iDepth,iTime] lt 0.27)then icolor=255
  endif else begin
   icolor = 255
  endelse

  ; plot data
  x0 = djulian - 1.d/48.d
  x1 = djulian
  y0 = soilHeight[iDepth]
  y1 = soilHeight[iDepth+1]
  polyfill, [x0,x1,x1,x0], [y0,y0,y1,y1], color=icolor

 endfor  ; soil layers
endfor  ; time

end
