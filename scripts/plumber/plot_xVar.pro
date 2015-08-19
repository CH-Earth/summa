pro plot_xVar

; define plotting parameters
window, 0, xs=2000, ys=1400, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3
!P.COLOR=0
erase, color=255
!P.MULTI=[0,5,4,0,0]

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
               ;'CABLE_2.0_SLI.vxh599_r553',  $
               'CHTESSEL',                   $
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
               'SUMMA.1.0']

; define colors
iColor=[80,140,210,250]

; define the number of models and sites
nModels = n_elements(model_names)
nSites  = n_elements(site_names)

; define file paths
file_path = '/d1/mclark/PLUMBER_data/model_output/'

; loop through sites
for iSite=0,nSites-1 do begin

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

  ; define the simulation time
  time0 = djulian_mod[0]
  time1 = djulian_mod[ntime_mod-1]

  ; make a base plot
  if(iModel eq 0)then begin
   plot, indgen(5), xrange=[time0,time1], yrange=[vmin,vmax], xstyle=9, ystyle=1, xticklen=(-0.02),$
    xtickformat=['label_date'], xticks=3, ytitle = cVarName, title=site_names[iSite]+'!C('+site_pfts[iSite]+')', $
    ymargin=[4,4], /nodata
   plots, [time0,time1], [vmax,vmax]
  endif

  ; plot the data
  oplot, djulian_mod, xVar, color=icolor[iModel]

  ; plot a legend
  if(iSite eq 0)then begin
   x0 = 0.05*(time1 - time0) + time0
   x1 = 0.40*(time1 - time0) + time0
   y0 = 0.70*(vmax - vmin) + vmin + 0.1*(vmax - vmin)*iModel
   y1 = y0 - 0.025*(vmax - vmin)
   plots, [x0,x1], [y0,y0], color=icolor[iModel]
   xyouts, x1 + 0.05*(time1 - time0), y1, model_names[iModel], charsize=2
  endif

 endfor  ; looping through models

endfor  ; looping through sites

; write figure
write_png, 'figures/plumber_'+cVarName+'.png', tvrd(true=1)

stop
end
