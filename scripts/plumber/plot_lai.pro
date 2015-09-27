pro plot_lai

; define plotting parameters
window, 0, xs=2000, ys=1400, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3
!P.COLOR=0
erase, color=255
!P.MULTI=[0,4,5,0,1]

; define the date format
dummy = label_date(date_format=['%D-%M!C%Y'])

; define the HRU
iHRU=0

; define variable
cVarName='LAI'

; define variable range
vmin=0
vmax=8

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

; define the model names
model_names = ['CABLE.2.0',                  $
               'SUMMA.1.0.exp.02.030',       $
               ;'CABLE_2.0_SLI.vxh599_r553',  $
               'CHTESSEL',                   $
               'SUMMA.1.0.exp.02.031',       $
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
               'SUMMA.1.0.exp.02.032',       $
               'SUMMA.1.0.exp.02.000']

; define colors
iColor=[40,80,120,160,210,250]
;iColor=[40,120,210]

; define the number of models and sites
nModels = n_elements(model_names)
nSites  = n_elements(site_names)

; define file paths
file_path = '/d1/mclark/PLUMBER_data/model_output/'

; loop through sites
for iSite=0,nSites-1 do begin

  ; define the tick names
  xtick_labels = [' ', strtrim(indgen(12)+1,2), ' ']

  ; make a base plot
  plot, indgen(5), xrange=[0,13], yrange=[vmin,vmax], xstyle=1, ystyle=1, xticklen=(-0.02),$
   xticks=13, ytitle = cVarName, title=site_names[iSite]+'!C('+site_pfts[iSite]+')', $
   ymargin=[4,4], xtickname=xtick_labels, /nodata

 ; get the filename for the Noah LAI
 ;lai_file = '/home/mclark/summa/settings/plumber/LAI_tables/time_parms_' + strtrim(iSite+1,2) + '.txt'

 ; define vectors for LAI and SAI
 ;xLAI = fltarr(12)
 ;xSAI = fltarr(12)

 ; read the LAI for the Noah simulations
 ;openr, lai_unit, lai_file, /get_lun
 ; readf, lai_unit, xLAI
 ; readf, lai_unit, xSAI
 ;free_lun, lai_unit
 
 ; plot it up
 ;oplot, indgen(12)+1, xLAI+xSAI, color=40, thick=3

 ; loop through the desired PLUMBER models
 for iModel=0,nModels-1 do begin

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
   xVar=reform(xVar)

  ; close the netcdf file
  ncdf_close, ncFileID

  ; get the time of the year
  caldat, djulian_mod, im, id, iyyy

  ; average the lai
  xLAI = fltarr(12)
  for imonth=1,12 do begin
   iMatch = where(im eq imonth, nMatch)
   if(nMatch gt 0)then xLAI[imonth-1] = mean(xVar[iMatch])
  endfor

  ; plot the data
  oplot, indgen(12)+1, xLAI, color=icolor[iModel]

  ; plot a legend
  if(iSite eq 0)then begin
   x0 = 0.5
   x1 = 4.0
   y0 = 0.70*(vmax - vmin) + vmin + 0.1*(vmax - vmin)*iModel
   y1 = y0 - 0.025*(vmax - vmin)
   plots, [x0,x1], [y0,y0], color=icolor[iModel]
   xyouts, x1 + 0.55, y1, model_names[iModel], charsize=2
  endif

 endfor  ; looping through models

 stop

endfor  ; looping through sites

; write figure
write_png, 'figures/plumber_'+cVarName+'-test.png', tvrd(true=1)

stop
end
