pro plot_paper2figure05

; define plotting parameters
window, 1, xs=1400, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=1.5
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,2,0,0]

; define the path and name of the graphics file
gpath = '/home/mclark/test_summa/summa/verification/zFigures/'
gname = 'Clark_et_al__WRR2015b_figure05.png'

; define the path of the validation data
valPath = '/home/mclark/test_summa/summa/testCases_data/validationData/'

; define list of obs files
file_list = [valPath+'Storck_9697_cutTreeData.txt', $
             valPath+'Storck_9798_cutTreeData.txt']

; define file path for model output
file_path = '/home/mclark/test_summa/summa/output/wrrPaperTestCases/figure05/'

; define file prefix for different experiments
file_pref = ['storckSite_spinup_hedpom', $
             'storckSite_spinup_storck']

; define file suffix
file_suff = ['9697', $
             '9798']

; define the number of tick marks
nxticks = [9,8]

; define initial HRU
iHRU = 0

; loop through files
for ifile=0,n_elements(file_suff)-1 do begin
 
 ; loop through experiments
 for iPrefix=0,n_elements(file_pref)-1 do begin

  ; *****
  ; (1) GET BASIC INFO FROM MODEL FILE...
  ; *************************************

  ; define file name
  filenm_mod = file_path + file_pref[iPrefix] + file_suff[ifile] + '.nc'
  print, filenm_mod

  ; open file
  nc_file = ncdf_open(filenm_mod, /nowrite)

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
   djulian_mod = bjulian + atime*aoff

   ; get the number of time elements
   ntime = n_elements(djulian_mod)-1

  ; close the NetCDF file
  ncdf_close, nc_file

  ; *****
  ; (2) PLOT OBSERVATIONS OF CANOPY STORAGE...
  ; ******************************************

  ; define filename
  filenm_obs = file_list[ifile]

  ; identify number of lines
  nlines = file_lines(filenm_obs)
 
  ; define data
  dataArr = dblarr(9,nlines)

  ; read data
  openr, in_unit, filenm_obs, /get_lun
   readf, in_unit, dataArr
  free_lun, in_unit

  ; get time
  iy  = long(dataArr[0,*])
  im  = long(dataArr[1,*])
  id  = long(dataArr[2,*])
  ih  = long(dataArr[3,*])
  imi = long(dataArr[4,*])

  ; get data
  belowCanopySWE = reform(dataArr[5,*])
  tree01_storage = reform(dataArr[6,*])
  tree02_storage = reform(dataArr[7,*])
  tree03_storage = reform(dataArr[8,*])

  ; get julian day
  djulian = julday(im,id,iy,ih,imi)

  ; define plot range
  if(ifile eq 0)then begin
   i_beg = 248
   i_end = 6728
  endif else begin
   i_beg = 624
   i_end = 6432
  endelse

  ; define the date format
  dummy = label_date(date_format=['%D %M!C%Y'])

  ; define x margin
  if(iprefix eq 0)then xmar=[8,2]
  if(iprefix eq 1)then xmar=[6,4]

  ; define ymargin
  if(ifile eq 0)then ymar=[4,2]
  if(ifile eq 1)then ymar=[6,0]

  ; make a base plot
  plot, djulian[i_beg:i_end], xrange=[djulian[i_beg],djulian[i_end]], yrange=[0,40], xstyle=1, ystyle=1, $
   xtickformat=['label_date'], xticks=nxticks[ifile], ytitle = 'Canopy interception (mm)', xmargin=xmar, ymargin=ymar, $
   /nodata

  ; define the color grey
  tvlct, r, g, b, /get
  r[1] = 180
  g[1] = 180
  b[1] = 180
  tvlct, r, g, b
 
  ; overplot data
  for itime=i_beg,i_end do begin
   if(tree01_storage[itime] gt 0.d)then $
    oplot, [djulian[itime],djulian[itime]], [0.d,tree01_storage[itime]], color=1
  endfor
 
  ; identify area of no data
  if(ifile eq 0)then ktime = 2250
  if(ifile eq 1)then ktime = 1265
  oplot, [djulian[ktime],djulian[ktime]], [0.d,40.d], color=0
  xyouts, mean([djulian[i_beg],djulian[ktime]]), 35.d, 'no data', alignment=0.5


 
  ; *****
  ; (3) PLOT MODEL SIMULATIONS...
  ; *****************************

  ; define colors
  icolor=[80,254]

  ; loop through HRUs
  for iHRU=0,1 do begin
 
   ; open netcdf file
   nc_file = ncdf_open(filenm_mod, /nowrite)
 
    ; get intercepted liquid water
    ivar_id = ncdf_varid(nc_file,'scalarCanopyLiq')
    ncdf_varget, nc_file, ivar_id, scalarCanopyLiq, offset=[iHRU,0], count=[1,ntime]
    scalarCanopyLiq = reform(scalarCanopyLiq)
 
    ; get intercepted ice
    ivar_id = ncdf_varid(nc_file,'scalarCanopyIce')
    ncdf_varget, nc_file, ivar_id, scalarCanopyIce, offset=[iHRU,0], count=[1,ntime]
    scalarCanopyIce = reform(scalarCanopyIce)
 
   ; close netcdf file
   ncdf_close, nc_file
 
   oplot, djulian_mod, scalarCanopyLiq+scalarCanopyIce, color=icolor[iHRU], thick=1.5

  endfor  ; (looping thru HRUs)

  ; define position of the legend
  y_top = 36.0
  y_gap =  2.0
  xline = [djulian[i_beg]+0.45d*(djulian[i_end]-djulian[i_beg]), $
           djulian[i_beg]+0.60d*(djulian[i_end]-djulian[i_beg])]
  ; define title for the legend
  jcolor=[80,254]
  if(iprefix eq 0)then xyouts, xline[0], y_top+1.5d, 'Hedstrom and Pomeroy (1998)', charsize=1.5
  if(iprefix eq 1)then xyouts, xline[0], y_top+1.5d, 'Andreadis et al. (2009)', charsize=1.5
  ; define text for the legend
  if(iprefix eq 0)then xtext = ['Branch capacity (cold) = 6.6 mm', 'Branch capacity (cold) = 13.2 mm']
  if(iprefix eq 1)then xtext = ['Branch capacity (cold) = 5 mm',   'Branch capacity (cold) = 10 mm']
  ; define the legend
  for ivar=0,1 do begin
   yline = y_top - float(ivar)*y_gap
   plots, xline, [yline,yline], color=jcolor[ivar]
   xyouts, xline[1]+1., yline-0.25, xtext[ivar], charsize=1.25
  endfor
  ; define the observations
  x1 = xline[0]
  x2 = xline[1]
  y1 = y_top - float(ivar)*y_gap - y_gap*0.65
  y2 = y_top - float(ivar)*y_gap + y_gap*0.35
  polyfill, [x1,x2,x2,x1], [y1,y1,y2,y2], color=1
  xyouts, x2+1., y1+0.75, 'Observed canopy SWE', charsize=1.25

 endfor  ; (looping through decisions)

endfor  ; looping through years

; make a figure
write_png, gpath+gname, tvrd(true=1)


stop
end
