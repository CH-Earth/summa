pro processPlumberInputData

; define the HRU index
ixHRU=1001

; define indices for the plant functional types
; NOTE: use USGS tables
ixCropland            = 3  ; Irrigated Cropland and Pasture
ixGrassland           = 7  ; Grassland
ixSavanna             = 10 ; Savanna
ixWoodySavanna        = 6  ; Cropland/Woodland Mosaic
ixDeciduousBroadleaf  = 11 ; Deciduous Broadleaf Forest
ixDeciduousNeedleleaf = 12 ; Deciduous Needleleaf Forest
ixEvergreenBroadleaf  = 13 ; Evergreen Broadleaf Forest
ixEvergreenNeedleleaf = 14 ; Evergreen Needleleaf Forest
ixMixedForest         = 15 ; Mixed Forest
ixPermanentWetland    = 17 ; Herbaceous Wetland

; define the path to the input files
inputData_path = '/home/mclark/summa/input/plumber/'

; define the settings path
settings_path = '/home/mclark/summa/settings/plumber/'

; define the path to the plumber data
plumber_path = '/d1/mclark/PLUMBER_data/'

; define the path to the site data
site_path =  plumber_path + 'site_data/met/'

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

; define plant functional types
site_pfts = [ixGrassland,          $       ; 'Amplero',     
             ixEvergreenNeedleleaf,$       ; 'Blodgett',    
             ixGrassland,          $       ; 'Bugac',       
             ixEvergreenNeedleleaf,$       ; 'ElSaler',    
             ixCropland,           $       ; 'ElSaler2',     
             ixEvergreenBroadleaf, $       ; 'Espirra',     
             ixGrassland,          $       ; 'FortPeck',    
             ixDeciduousBroadleaf, $       ; 'Harvard',     
             ixDeciduousBroadleaf, $       ; 'Hesse',       
             ixWoodySavanna,       $       ; 'Howard',      
             ixEvergreenNeedleleaf,$       ; 'Howlandm',    
             ixEvergreenNeedleleaf,$       ; 'Hyytiala',    
             ixSavanna,            $       ; 'Kruger',      
             ixEvergreenNeedleleaf,$       ; 'Loobos',      
             ixPermanentWetland,   $       ; 'Merbleue',    
             ixWoodySavanna,       $       ; 'Mopane',      
             ixEvergreenBroadleaf, $       ; 'Palang',      
             ixMixedForest,        $       ; 'Sylvania',    
             ixEvergreenBroadleaf, $       ; 'Tumba',       
             ixDeciduousBroadleaf  ]       ; 'UniMich'      

; define the site attributes
site_attr = ['latitude',    $
             'longitude',   $
             'elevation',   $
             'reference_height']

; define the forcing variables
site_vars = ['Rainf',    $
             'SWdown',   $
             'LWdown',   $
             'Tair',     $
             'Wind',     $
             'PSurf',    $
             'Qair'      ]

; define the variable types
ix3d = 0
ix4d = 1

; define the type of forcing variables
var_type = [ix3d,  $  ; 'Rainf',    
            ix3d,  $  ; 'SWdown',   
            ix3d,  $  ; 'LWdown',   
            ix4d,  $  ; 'Tair',     
            ix4d,  $  ; 'Wind',     
            ix3d,  $  ; 'PSurf',    
            ix4d   ]  ; 'Qair'      

; define the number of sites, attributes and variables
nSites  = n_elements(site_names)
nAtts   = n_elements(site_attr)
nVars   = n_elements(site_vars)

; define the vector of attributes and variables
xAttributes = dblarr(nAtts)
xForcing    = dblarr(nVars)

; define un-used attributes
HRUarea        = 1.d
tan_slope      = 0.2d
contourLength  = 100.d
soilTypeIndex  = 8
slopeTypeIndex = 1
downHRUindex   = 0

; define the header for the local attributed file
cHead_attr = 'hruIndex    HRUarea   latitude  longitude  elevation  tan_slope  contourLength' + $
             '  mHeight   vegTypeIndex  soilTypeIndex slopeTypeIndex   downHRUindex'

; loop through sites
for iSite=0,nSites-1 do begin

 ; print progress
 print, iSite+1, ': ', site_names[iSite]

 ; make a directory for a given site
 ;if(file_test(settings_path + site_names[iSite]) eq 0)then $
 ; spawn, 'mkdir ' + settings_path + site_names[iSite]

 ; define filename for the input data
 inputData_filename = inputData_path + 'inputData_' + site_names[iSite] + '.txt'

 ; define filename for the attributes
 attr_filename1 = settings_path + 'attributes/summa_zLocalAttributes_' + site_names[iSite] + '.txt'
 attr_filename2 = settings_path + 'attributes_siteSpecific/summa_zLocalAttributes_' + site_names[iSite] + '.txt'

 ; open up the input data file
 openw, outunit_data, inputData_filename, /get_lun

 ; open up the attributes files
 openw, outunit_att1, attr_filename1, /get_lun
 openw, outunit_att2, attr_filename2, /get_lun

 ; write the header for the local attributes file
 printf, outunit_att1, cHead_attr
 printf, outunit_att2, cHead_attr

 ; define filename
 site_filename = site_path + site_names[iSite] + 'Fluxnet.1.4_met.nc'

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
  ntime = n_elements(djulian)-1

  ; read in the site attributes
  for iAtt=0,nAtts-1 do begin
   ivar_id = ncdf_varid(nc_file,site_attr[iAtt])
   ncdf_varget, nc_file, ivar_id, xAtt
   xAttributes[iAtt] = reform(xAtt)
  endfor  ; attributes

  ; print the attributes using the USGS tables
  printf, outunit_att1, ixHRU, HRUarea, xAttributes[0:2], tan_slope, contourLength, xAttributes[3], $
           site_pfts[iSite], soilTypeIndex, slopeTypeIndex, downHRUindex, format='(i8,1x,4(f10.3,1x),f10.5,1x,f14.5,1x,f8.2,1x,4(i14,1x))'

  ; print the attributes using the site-specific tables
  printf, outunit_att2, ixHRU, HRUarea, xAttributes[0:2], tan_slope, contourLength, xAttributes[3], $
           iSite+1, soilTypeIndex, slopeTypeIndex, downHRUindex, format='(i8,1x,4(f10.3,1x),f10.5,1x,f14.5,1x,f8.2,1x,4(i14,1x))'

  ; loop through time
  for itime=0,ntime-1 do begin

   ; get the dates
   caldat, djulian[itime], im, id, iyyy, ih, imi, asec

   ; get the variables
   for iVar=0,nVars-1 do begin
    ivar_id = ncdf_varid(nc_file,site_vars[iVar])
    if(var_type[iVar] eq ix3d)then ncdf_varget, nc_file, ivar_id, xVar, offset=[0,0,itime], count=[1,1,1]
    if(var_type[iVar] eq ix4d)then ncdf_varget, nc_file, ivar_id, xVar, offset=[0,0,0,itime], count=[1,1,1,1]
    xForcing[iVar] = reform(xVar)
   endfor

   ; print data
   printf, outunit_data, iyyy, im, id, ih, imi, asec, xForcing, format='(i4,1x,4(i2,1x),f6.1,1x,e14.4,1x,5(f10.3,1x),e12.3)'

  endfor  ; looping through time

 ; close netcdf file
 ncdf_close, nc_file

 ; free up logical units for the output files
 free_lun, outunit_att1
 free_lun, outunit_att2
 free_lun, outunit_data

endfor  ; looping through the sites


stop
end
