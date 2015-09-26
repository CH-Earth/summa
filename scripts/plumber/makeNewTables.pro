pro makeNewTables

; define the HRU
iHRU=0

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

; define look-up table to match soutthern hemisphere sites to Northern hemisphere
ixMatchSH = [7,8,9,10,11,12,1,2,3,4,5,6]
     
; define the model names
model_names = ['CABLE.2.0','CHTESSEL','SUMMA.1.0.exp.02.000']

; define names for the tables
table_names = ['plumberCABLE','plumberCHTESSEL','plumberSUMMA']

; define name for latitude
lat_name = ['latitude','lat','latitude']

; define the number of models and sites
nModels = n_elements(model_names)
nSites  = n_elements(site_names)

; define file paths
model_path = '/d1/mclark/PLUMBER_data/model_output/'

; define table
table_path = '/home/mclark/summa/settings/plumber/'

; loop through the desired PLUMBER models
for iModel=0,nModels-1 do begin

 ; *****
 ; * GET THE LAI...
 ; ****************

 ; define the latitude
 xLatitude = fltarr(nSites)

 ; define the SAI, LAI, and VAI
 xSAI = fltarr(12,nSites) ; stem area index
 xLAI = fltarr(12,nSites) ; leaf area index
 xVAI = fltarr(12,nSites) ; vegetation area index (LAI+SAI)

 ; loop through the sites
 for iSite=0,nSites-1 do begin

  ; define the file name
  file_name = model_names[iModel] + '/' + model_names[iModel] + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'

  ; open files
  ncFileID = ncdf_open(model_path+file_name, /nowrite)

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

   ; get the LAI (actually VAI)
   ivar_id = ncdf_varid(ncFileID,'LAI')
   ncdf_varget, ncFileID, ivar_id, xVar
   xVar=reform(xVar)

   ; get the latitude
   ivar_id = ncdf_varid(ncFileID,lat_name[iModel])
   ncdf_varget, ncFileID, ivar_id, xLat
   xLat=reform(xLat)

  ; close the netcdf file
  ncdf_close, ncFileID

  ; get the time of the year
  caldat, djulian_mod, im, id, iyyy

  ; average the lai
  for imonth=1,12 do begin
   iMatch = where(im eq imonth, nMatch)
   if(nMatch gt 0)then xVAI[imonth-1,iSite] = mean(xVar[iMatch])
  endfor

  ; save the latitude
  xLatitude[iSite] = xLat

 endfor  ; looping through sites

 ; *****
 ; * READ/WRITE NOAH-MP TABLE...
 ; *****************************

 ; define the name of the NoahMP table
 tableNoahMP = table_path + 'MPTABLE.TBL'

 ; define table type
 table_type = 'USGS'

 ; define filename for a new table
 tableNoahMP_new = table_path + table_names[iModel] + '_MPTABLE.TBL'

 ; open file for writing
 openw, tableUnitOut, tableNoahMP_new, /get_lun

 ; write header
 printf, tableUnitOut, ' '
 printf, tableUnitOut, ' '
 printf, tableUnitOut, '&noah_mp_' + table_names[iModel] + '_veg_categories'
 printf, tableUnitOut, ' VEG_DATASET_DESCRIPTION = "LAI used in the ' + model_names[iModel] + ' PLUMBER simulations"'
 printf, tableUnitOut, ' NVEG = ' + strtrim(nSites,2)
 printf, tableUnitOut, '/'
 printf, tableUnitOut, '&noah_mp_' + table_names[iModel] + '_parameters'
 printf, tableUnitOut, ' ! NVEG = ' + strtrim(nSites,2)

 ; print the site names
 for iSite=0,nSites-1 do begin
  printf, tableUnitOut, ' ! ' + string(iSite+1, format='(i2)') + ': ' + site_names[iSite]
 endfor  ; looping through sites
 printf, tableUnitOut, ' '

 ; define a character string
 cLine=''

 ; intialize flag for finding the table
 foundTable=0
 startTable=0
 doneHeader=0

 ; define counter
 iCount=0

 ; get the parameter vector for each site
 yVector = strarr(nSites)

 ; open table for reading
 openr, tableUnitIn, tableNoahMP, /get_lun

  ; read to the end
  while ~eof(tableUnitIn) do begin

   ; read a line of data
   readf, tableUnitIn, cLine

   ; identify when we have found the data
   if(strmatch(cLine,'*"'+table_type+'"*') eq 1)then foundTable=1
   if(foundTable eq 0)then continue
   
   ; get the start of the table
   if(strmid(cLine,0,1) eq '/')then startTable=1
   if(startTable eq 0)then continue

   ; exit if get the end of the table
   if(iCount gt 0 and strmid(cLine,0,1) eq '/')then break

   ; ignore comments and blank lines
   if(strmid(strtrim(cLine,2),0,1) eq '!' or strlen(cLine) eq 0)then continue

   ; ignore tags
   if(strmid(strtrim(cLine,2),0,1) eq '&' or strmid(strtrim(cLine,2),0,1) eq '/')then continue

   ; deal with scalars (no commas)
   if(strmatch(cLine,'*,*') eq 0)then begin
    cTemp = strsplit(cLine, '=', /extract)
    printf, tableUnitOut, cTemp[0] + ' = -999'

   ; deal with vectors
   endif else begin

    ; print the header of the table
    if(doneHeader eq 0)then begin
     printf, tableUnitOut, ' '
     printf, tableUnitOut, ' !-------------------------------------------------------------------------------------------------------------------------------------------------'
     printf, tableUnitOut, ' !     ' + strjoin(string(indgen(nsites)+1,format='(i7)'))
     printf, tableUnitOut, ' !-------------------------------------------------------------------------------------------------------------------------------------------------'
     doneHeader=1
    endif

    ; extract the name
    if(strmatch(cLine,'*=*') eq 1)then begin
     cTemp = strsplit(cLine, '=', /extract)
     cName = cTemp[0]
     cData = cTemp[1]
     cLabel = cName + '='
     cNameSave = cName
    endif else begin
     cName = ' '
     cData = cLine
     cLabel = '        '
    endelse

    ; don't write LAI
    if(strtrim(cNameSave,2) eq 'LAIM')then continue

    ; add some comments
    if(strtrim(cName,2) eq 'RHOL' or strtrim(cName,2) eq 'RHOS' or strtrim(cName,2) eq 'TAUL' or strtrim(cName,2) eq 'TAUS')then begin
     printf, tableUnitOut, ' '
     printf, tableUnitOut, ' ! Row 1:  Vis'
     printf, tableUnitOut, ' ! Row 2:  Near IR'
    endif

    ; add some more comments
    if(strtrim(cName,2) eq 'EPS')then begin
     printf, tableUnitOut, ' '
     printf, tableUnitOut, ' ! Five types, one row for each type.'
    endif
 
    ; add a blank line
    if(strtrim(cName,2) eq 'XL' or strtrim(cName,2) eq 'LTOVRC' or strtrim(cName,2) eq 'SAIM')then printf, tableUnitOut, ' '

    ; extract the vector
    xVector = strsplit(cData,',',/extract)

    ; get the vector for each site
    for iSite=0,nSites-1 do yVector[iSite] = xVector[site_pfts[iSite]-1]

    ; make it into a string
    cData = strjoin(yVector,',', /single)

    ; save SAI
    if(strtrim(cNameSave,2) eq 'SAIM')then begin
     if(strtrim(cName,2) eq 'SAIM')then iMonth=0 else iMonth=iMonth+1
     xSAI[iMonth,*] = float(strsplit(cData, ',', /extract))
    endif

    ; print the LAI
    if(strtrim(cName,2) eq 'SLAREA')then begin

     ; compute the LAI
     for iSite=0,nSites-1 do begin
      for iMonth=0,11 do begin

       ; define the month index (different for Southern Hemisphere sites
       if(xLatitude[iSite] lt 0.d)then jmonth=ixMatchSH[iMonth]-1 else jmonth=imonth

       ; subtract the stem area index
       ; NOTE: use of jmonth for xVAI
       xLAI[iMonth,iSite] = max([0.d, xVAI[jMonth,iSite] - xSAI[iMonth,iSite]])

      endfor
     endfor

     ; print the LAI
     printf, tableUnitOut, ' '
     for iMonth=0,11 do begin
      if(iMonth eq 0)then cID=' LAIM  =' else cID='        '
      printf, tableUnitOut, cID + strjoin(string(xLAI[iMonth,*],format='(f6.2)'),',',/single) + ','
     endfor
     printf, tableUnitOut, ' '

    endif ; if printing the LAI

    ; print the data
    printf, tableUnitOut, cLabel + cData + ','

   endelse  ; reading vectors

   ; increment counter
   iCount = iCount+1

  endwhile  ; looping through file

  ; define the end of the table
  printf, tableUnitOut, '/'

 ; free up file unit
 free_lun, tableUnitIn
 free_lun, tableUnitOut

 ; list the file
 spawn, 'cat ' + tableNoahMP_new

 ; *****
 ; * READ/WRITE VEG TABLE...
 ; *****************************

 ; define the name of the NoahMP table
 tableVeg = table_path + 'VEGPARM.TBL'

 ; define table type
 table_type = 'USGS'

 ; define filename for a new table
 tableVeg_new = table_path + table_names[iModel] + '_VEGPARM.TBL'

 ; intialize flag for finding the table
 foundTable=0
 startTable=0
 doneHeader=0

 ; define counter
 iCount=0

 ; open file for writing
 openw, tableUnitOut, tableVeg_new, /get_lun

 ; write header
 printf, tableUnitOut, 'Vegetation Parameters'
 printf, tableUnitOut, table_names[iModel]

 ; open table for reading
 openr, tableUnitIn, tableVeg, /get_lun

 ; read to the end
 while ~eof(tableUnitIn) do begin

  ; read line
  readf, tableUnitIn, cLine

  ; start once we have found the table
  if(strtrim(cLine,2) eq table_type)then begin

   ; read, modify and write the header
   readf, tableUnitIn, cLine               ; read header
   cTemp    = strsplit(cLine,',',/extract) ; split header
   nTypes   = fix(cTemp[0])
   cTemp[0] = strtrim(nSites,2)            ; modify the number of sites
   cHead    = strjoin(cTemp,',', /single)  ; recreate the header
   printf, tableUnitOut, cHead             ; print the header

   ; get an array of veg Types
   vegTypes = strarr(nTypes)

   ; read the veg types
   for iType=0,nTypes-1 do begin
    readf, tableUnitIn, cLine
    if(iType+1 lt 10)then begin ; just fixing the format
     cData    = strsplit(cLine, ',', /extract)
     cData[1] = strmid(cData[1],1)
     cLine    = strjoin(cData,',', /single)
    endif
    vegTypes[iType] = cLine
   endfor

   ; write the veg types
   for iSite=0,nSites-1 do begin
    cLine    = vegTypes[site_pfts[iSite]-1]      ; get the veg data for a given site
    cTemp    = strsplit(cLine, "'", /extract)    ; split out the name of the veg type
    cData    = strsplit(cTemp[0], ',', /extract) ; split out the data
    cData[0] = strtrim(iSite+1,2)                ; overwrite the site index
    if(iSite+1 lt 10)then cData[1]=' '+cData[1]  ; fix the format
    cLine    = strjoin(cData,',', /single)       ; recreate the data line
    printf, tableUnitOut, cLine + "     '" + site_names[iSite] + "'"
   endfor

   ; set flag to found
   foundTable=1
   continue

  endif  ; if found the table

  ; keep trying until we find the table
  if(foundTable eq 0)then continue

  ; break if we get to the next table
  if(foundTable eq 1 and strtrim(cLine,2) eq 'Vegetation Parameters')then break

  ; default: print out the line
  printf, tableUnitOut, cLine

 endwhile  ; looping through input table

 ; free up file unit
 free_lun, tableUnitIn
 free_lun, tableUnitOut

 ; list the file
 spawn, 'cat ' + tableVeg_new

endfor  ; looping through models


stop
end
