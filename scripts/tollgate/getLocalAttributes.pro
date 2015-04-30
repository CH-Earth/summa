pro getLocalAttributes

; define plotting parameters
window, 0, xs=1500, ys=750, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2.5
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,1]

; define x range
xmin = -1301000.d
xmax = -1285000.d

; define y range
ymin = 190000.d
ymax = 206000.d

; define x range
x1 = xmin + 0.00*(xmax - xmin)
x2 = xmin + 1.00*(xmax - xmin)

; define y range
y1 = ymin + 0.00*(ymax - ymin)
y2 = ymin + 1.00*(ymax - ymin)

; *****
; (1) GET THE INFORMATION FROM THE NETWORK TOPOLOGY FILE...
; *********************************************************

; declare network topology file
filenm = '/home/mclark/summa/ancillary_data/tollgate/Reynolds_Network_Topology.nc'
print, filenm

; open file
nc_file = ncdf_open(filenm, /nowrite)

 ; get the hruid
 ivar_id = ncdf_varid(nc_file,'hruid')
 ncdf_varget, nc_file, ivar_id, ixHRU

 ; get the basin area
 ivar_id = ncdf_varid(nc_file,'Basin_Area')
 ncdf_varget, nc_file, ivar_id, Basin_Area

 ; get the upstream area
 ivar_id = ncdf_varid(nc_file,'UpArea')
 ncdf_varget, nc_file, ivar_id, UpArea

 ; get the average elevation
 ivar_id = ncdf_varid(nc_file,'Elev_Avg')
 ncdf_varget, nc_file, ivar_id, Elev_Avg

 ; get the elevation at the downstream end of the stream segment
 ivar_id = ncdf_varid(nc_file,'BotElev')
 ncdf_varget, nc_file, ivar_id, BotElev

 ; get the index of the downstream segment
 ivar_id = ncdf_varid(nc_file,'tosegment')
 ncdf_varget, nc_file, ivar_id, ixDownstream

; close file
ncdf_close, nc_file

; define the desired reach
ixDesire = 281  ; the reach that includes the tollgate gauge

; define the total area above tollgate
totalArea = UpArea[ixDesire]

; define the number of reaches
nRch = n_elements(UpArea)

; define the reaches that flow through the tollgate reach
ixTollgate = intarr(nRch)

; identify if a given reach is desired
for iRch=0,nRch-1 do begin

 ; initialize the first downstream reach
 ixTrial = ixDownstream[iRch]

 ; climb down the network and check if we pass through ixDesire
 while(ixTrial ne -1) do begin
  if(ixTrial eq ixDesire+1)then ixTollgate[iRch]=1
  jxTrial = where(ixTrial eq ixHRU, nMatch)
  if(nMatch ne 1)then stop, 'no match'
  ixTrial = ixDownstream[jxTrial[0]]
 endwhile

endfor  ; looping through reaches

; add ixDesire
ixTollgate[ixDesire]=1

; get the number of basins
nBasins = total(ixTollgate)

print, 'botElev   = ', botElev[ixDesire]
print, 'totalArea = ', totalArea + Basin_Area[ixDesire]
print, 'nBasins   = ', nBasins

; *****
; (2) READ IN THE ASCII GRIDS FOR VEG...
; **************************************

; define the header
cHead=''
nHead=6

; define the file
filenm = '/home/mclark/summa/ancillary_data/tollgate/tollgate_veg_30m_dat_lcc.txt'

; open file for reading
openr, in_unit, filenm, /get_lun

 ; loop through header lines
 for iHead=0,nHead-1 do begin
  ; read header
  readf, in_unit, cHead
  cData = strsplit(cHead,' ',/extract)
  ; extract grid info
  case cData[0] of
   'ncols':        nCols     = long(cData[1])
   'nrows':        nRows     = long(cData[1])
   'xllcorner':    xll       = float(cData[1])
   'yllcorner':    yll       = float(cData[1])
   'cellsize':     cSize     = float(cData[1])
   'NODATA_value': ixMissing = long(cData[1]) 
   else: stop, 'unable to find header value'
  endcase
 endfor  ; end looping through header

 ; extract grid
 ixVeg = intarr(nCols,nRows)
 readf, in_unit, ixVeg
 ixVeg = temporary(reverse(ixVeg,2))

; close file
free_lun, in_unit 

; define x and y coordinates
xVeg = xll + dindgen(nCols)*cSize
yVeg = yll + dindgen(nRows)*cSize

; make a base plot
plot, indgen(5), xrange=[x1,x2], yrange=[y1,y2], $
 xmargin=[1,1], ymargin=[1,1], xstyle=13, ystyle=13, /nodata

; plot image
for iCol=0,nCols-2 do begin
 for jRow=0,nRows-2 do begin
  ; check for valid data
  if(ixVeg[iCol,jRow] ge 0)then begin
   ; define x-coordinates
   xx1 = xVeg[iCol]
   xx2 = xVeg[iCol+1]
   ; define y-coordinates
   yy1 = yVeg[jRow]
   yy2 = yVeg[jRow+1]
   ; define the color
   case ixVeg[iCol,jRow] of
    1: icolor=70  ; riparian
    2: icolor=160 ; fir
    3: icolor=120 ; aspen
    4: icolor=180 ; sage
    5: icolor=190 ; low sage
    6: icolor=190 ; grass
    else: stop, 'unable to identify the veg type'
   endcase
   ; plot the pixel
   polyfill, [xx1,xx2,xx2,xx1], [yy1,yy1,yy2,yy2], color=icolor
  endif  ; if valid data
 endfor
endfor

; define valid pixels
iValid = where(ixVeg ge 0)

; check percent in each category
for iCat=0,10 do begin
 iMatch = where(ixVeg[iValid] eq iCat, nMatch)
 print, 'category ', iCat, nMatch, format='(a,1x,i2,1x,i10,1x)'
endfor


; *****
; (3) PLOT UP THE SUB-BASINS...
; *****************************

; Define the shapefile
file_name = '/home/mclark/summa/ancillary_data/tollgate/Catchment2_Project.shp'

; Open the Shapefile
myshape=OBJ_NEW('IDLffShape', file_name)

 ; Get the number of entities so we can parse through them
 myshape->IDLffShape::GetProperty, N_ENTITIES=num_ent

 ; Parsing through the entities
 FOR ix=0, (num_ent-1) DO BEGIN
  ; Get the Attributes for entity x
  attr = myshape->IDLffShape::GetAttributes(ix)
  ; Get the attribute info
  myshape->GetProperty, ATTRIBUTE_INFO=attr_info
  ; Get entity
  ent=myshape->IDLffShape::GetEntity(ix)
  ; save x and y
  x = reform((*ent.vertices)[0,*])
  y = reform((*ent.vertices)[1,*])
  ; get the number of vertices
  nVerts = ent[0].n_vertices
  ; get the number of parts
  nParts = ent[0].n_parts
  ; save the indices that define the parts of the shape
  iStart = *ent[0].parts
  if(nParts gt 1)then begin
   iCount = [iStart[1:nParts-1] - iStart[0:nParts-2], nVerts - iStart[nParts-1]]
  endif else begin
   iCount = [nVerts]
  endelse
  ; Clean-up of pointers
  myshape->IDLffShape::DestroyEntity, ent
  ; plot data
  if(ixTollgate[ix] eq 1)then begin
   ; plot data
   for iPart=0,nParts-1 do begin
    xx = x[iStart[iPart]:iStart[iPart]+iCount[iPart]-1]
    yy = y[iStart[iPart]:iStart[iPart]+iCount[iPart]-1]
    plots, xx, yy
   endfor
   ;if(ixTollgate[ix] eq 1)then plots, [x,x[0]], [y,y[0]], color=icol
  endif
 ENDFOR  ; parsing through the entities

; Close the Shapefile
OBJ_DESTROY, myshape

; *****
; (3) READ IN THE ASCII GRIDS FOR VEG2POLYGON CORRESPONDENCE...
; *************************************************************

; define the file
filenm = '/home/mclark/summa/ancillary_data/tollgate/veg2polygon_mapping_lcc.txt'

; open file for reading
openr, in_unit, filenm, /get_lun

 ; loop through header lines
 for iHead=0,nHead-1 do begin
  ; read header
  readf, in_unit, cHead
  cData = strsplit(cHead,' ',/extract)
  ; extract grid info
  case cData[0] of
   'ncols':        nCols     = long(cData[1])
   'nrows':        nRows     = long(cData[1])
   'xllcorner':    xll       = float(cData[1])
   'yllcorner':    yll       = float(cData[1])
   'cellsize':     cSize     = float(cData[1])
   'NODATA_value': ixMissing = long(cData[1])   
   else: stop, 'unable to find header value'
  endcase
 endfor  ; end looping through header

 ; extract grid
 ixVeg2Poly = intarr(nCols,nRows)
 readf, in_unit, ixVeg2Poly
 ixVeg2Poly = temporary(reverse(ixVeg2Poly,2))

; close file
free_lun, in_unit

; *****
; (4) PUT DATA INTO THE POLYGON ARRAY...
; **************************************

; define variables
xLat  =  43.2d
xLon  = 243.2d
slope =   0.2d
cLen  = 100.d

; define indices
ixSoil  = 8
ixSlope = 1
ixDowns = 0

; define array dimensioned by the number of polygons in the shape file
ixMode = intarr(nRch)

; put data into the polygon array
for iRch=0,nRch-1 do begin
 if(ixTollgate[iRch] eq 1)then begin

  ; get veg grid cells within each polygon
  iMatch = where(ixHRU[iRch] eq ixVeg2Poly, nMatch)
  if(nMatch eq 0)then stop, 'expect some grid cells in each polygon'

  ; get the mode
  ixHist = histogram(ixVeg[iMatch], locations=ixTable)
  ixMax  = max(ixHist,ixLoc)
  ixMode[iRch] = ixTable[ixLoc]

  ; get the fraction in each cell
  aFrac = double(ixHist)/double(nMatch)
  iGrab = where(aFrac gt 0.05d, nGrab)
  if(nGrab eq 0)then stop, 'expect some fraction greater than 5%'

  ; make the local attributes file
  file_path = '/home/mclark/summa/settings/tollgate/localAttributes/'
  file_name = file_path + 'snow_zLocalAttributes_basin_' + strtrim(string(ixHRU[iRch],format='(i4.4)'),2) + '.txt'
  openw, out_unit, file_name, /get_lun

  ; print header
  print, 'hruIndex    HRUarea   latitude  longitude  elevation  tan_slope  contourLength  mHeight  vegTypeIndex soilTypeIndex slopeTypeIndex downHRUindex', format='(a)'
  printf, out_unit, 'hruIndex    HRUarea   latitude  longitude  elevation  tan_slope  contourLength  mHeight  vegTypeIndex soilTypeIndex slopeTypeIndex downHRUindex', format='(a)'

  ; loop through HRUs
  for jGrab=0,nGrab-1 do begin

   ; get the HRU area
   xArea = Basin_Area[iRch] * aFrac[iGrab[jGrab]] / total(aFrac[iGrab])

   ; get the Seyfried veg type
   ixType = ixTable[iGrab[jGrab]]

   ; get the measurement height
   case ixType of
    1: xHeight =  3.d  ; riparian
    2: xHeight = 22.d  ; fir
    3: xHeight = 19.d  ; aspen
    4: xHeight =  3.d  ; sage
    5: xHeight =  3.d  ; low sage
    6: xHeight =  3.d  ; grass
    else: stop, 'unable to identify the veg type'  
   endcase

   ; get the Modis veg type
   case ixType of
    1: ixModis = 18 ; riparian
    2: ixModis = 13 ; fir
    3: ixModis = 11 ; aspen
    4: ixModis =  8 ; sage
    5: ixModis =  9 ; low sage
    6: ixModis =  9 ; grass
    else: stop, 'unable to identify the veg type'
   endcase

   ; write the attributes to the screen
   print, jGrab+1001, xArea, xLat, xLon, Elev_Avg[iRch], slope, cLen, xHeight, ixModis, ixSoil, ixSlope, ixDowns, $
    format='(i8,1x,4(f10.1,1x),f10.5,1x,f14.2,1x,f8.2,1x,2(i13,1x),i14,1x,i12)'

   ; write the attributes to the file
   printf, out_unit, jGrab+1001, xArea, xLat, xLon, Elev_Avg[iRch], slope, cLen, xHeight, ixModis, ixSoil, ixSlope, ixDowns, $
    format='(i8,1x,4(f10.1,1x),f10.5,1x,f14.2,1x,f8.2,1x,2(i13,1x),i14,1x,i12)'

  endfor  ; (looping through HRUs)

  ; free up the file unit
  free_lun, out_unit

  print, ixTable, format='(10(i4,1x))'
  print, ixHist, format='(10(i4,1x))'
  print, '**', ixMode[iRch], nGrab

 endif else begin
  ixMode[iRch] = -99
 endelse

endfor

; *****
; (5) PLOT UP THE SUB-BASINS...
; *****************************

; make a base plot
plot, indgen(5), xrange=[x1,x2], yrange=[y1,y2], $
 xmargin=[1,1], ymargin=[1,1], xstyle=13, ystyle=13, /nodata

; Define the shapefile
file_name = '/home/mclark/summa/ancillary_data/tollgate/Catchment2_Project.shp'

; Open the Shapefile
myshape=OBJ_NEW('IDLffShape', file_name)

 ; Get the number of entities so we can parse through them
 myshape->IDLffShape::GetProperty, N_ENTITIES=num_ent

 ; Parsing through the entities
 FOR ix=0, (num_ent-1) DO BEGIN
  ; Get the Attributes for entity x
  attr = myshape->IDLffShape::GetAttributes(ix)
  ; Get the attribute info
  myshape->GetProperty, ATTRIBUTE_INFO=attr_info
  ; Get entity
  ent=myshape->IDLffShape::GetEntity(ix)
  ; save x and y
  x = reform((*ent.vertices)[0,*])
  y = reform((*ent.vertices)[1,*])
  ; get the number of vertices
  nVerts = ent[0].n_vertices
  ; get the number of parts
  nParts = ent[0].n_parts
  ; save the indices that define the parts of the shape
  iStart = *ent[0].parts
  if(nParts gt 1)then begin
   iCount = [iStart[1:nParts-1] - iStart[0:nParts-2], nVerts - iStart[nParts-1]]
  endif else begin
   iCount = [nVerts]
  endelse
  ; Clean-up of pointers
  myshape->IDLffShape::DestroyEntity, ent
  ; plot data
  if(ixTollgate[ix] eq 1)then begin
   ; define the color
   case ixMode[ix] of
    1: icolor=70  ; riparian
    2: icolor=160 ; fir
    3: icolor=120 ; aspen
    4: icolor=180 ; sage
    5: icolor=190 ; low sage
    6: icolor=190 ; grass
    else: stop, 'unable to identify the veg type'
   endcase
   ; plot data
   for iPart=0,nParts-1 do begin
    xx = x[iStart[iPart]:iStart[iPart]+iCount[iPart]-1]
    yy = y[iStart[iPart]:iStart[iPart]+iCount[iPart]-1]
    polyfill, [xx,xx[0]], [yy,yy[0]], color=icolor
    plots, xx, yy
   endfor
   ;if(ixTollgate[ix] eq 1)then plots, [x,x[0]], [y,y[0]], color=icol
  endif
 ENDFOR  ; parsing through the entities

; Close the Shapefile
OBJ_DESTROY, myshape

; *****
; (4) PLOT UP THE DRAINAGE NETWORK...
; ***********************************

; Define the shapefile
file_name = '/home/mclark/summa/ancillary_data/tollgate/DrainageLine2_Project.shp'

; Open the Shapefile
myshape=OBJ_NEW('IDLffShape', file_name)

 ; Get the number of entities so we can parse through them
 myshape->IDLffShape::GetProperty, N_ENTITIES=num_ent

 ; Parsing through the entities
 FOR ix=0, (num_ent-1) DO BEGIN
  ; Get the Attributes for entity x
  attr = myshape->IDLffShape::GetAttributes(ix)
  ; Get the attribute info
  myshape->GetProperty, ATTRIBUTE_INFO=attr_info
  ; Get entity
  ent=myshape->IDLffShape::GetEntity(ix)
  ; save x and y
  x = reform((*ent.vertices)[0,*])
  y = reform((*ent.vertices)[1,*])
  ; Clean-up of pointers
  myshape->IDLffShape::DestroyEntity, ent
  ; get the desired segID
  iSegID = attr.ATTRIBUTE_6
  ; plot data
  if(iSegID le nRch)then begin
   ithick = 5.*(upArea[iSegID-1]/totalArea)
   if(ixTollgate[iSegID-1] eq 1)then plots, x, y, thick=ithick
  endif
 ENDFOR  ; parsing through the entities

; Close the Shapefile
OBJ_DESTROY, myshape


stop
end
