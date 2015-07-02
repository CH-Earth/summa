pro plot_paper2figure06

; define plotting parameters
window, 0, xs=1000, ys=500, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=1.5
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,1,0,0]

; define constants
Tfreeze = 273.16
iden_ice= 917.

nsites = 2

; define the path and name of the graphics file
gpath = '/home/mclark/test_summa/summa/verification/zFigures/'
gname = 'Clark_et_al__WRR2015b_figure06.png'

; define the path of the validation data
valPath = '/home/mclark/test_summa/summa/testCases_data/validationData/'

; define datafiles with model validation
valFile = [valPath+'ReynoldsCreek_valData.nc', $
           valPath+'senatorBeck_SASP_1hr.nc']

; define name of validation data
valName = ['zs_sheltered','snowDepth']

; define multiplier
valMult = [0.01d,1.d]  ; convert from cm to m, negative because coordinate variable is positive downwards with zero at the soil surface

; define file path
file_path = '/home/mclark/test_summa/summa/output/wrrPaperTestCases/figure06/'

; define file prefix
file_pref = 'albedoTest_'

; define site name
site_name = ['_reynolds', $
             '_senator']

; define the file suffix
cSuffix = ['VariableDecayRate','ConstantDecayRate']

; define colors for each experiment (defined by file suffix)
iColor=[80,250]

; define the plot title
ptitle = ['Reynolds Mountain East', 'Senator Beck']

; define the HRU
iHRU=0

; define year
iyear = [2005,2010]

; define x margin
xmar=[8,4]

; define ytitle
ytit='Snow depth (m)'

; loop through files
for ifile=0,nsites-1 do begin

 ; define water year
 cWaterYear = strtrim(iyear[ifile],2)+'-'+strtrim(iyear[ifile]+1,2)

 ; declare file
 filenm = file_path+file_pref+cWaterYear+site_name[ifile]+cSuffix[0]+'.nc'

 ; *****
 ; * GET BASIC INFO FROM MODEL FILE...
 ; ***********************************

 ; open file
 nc_file = ncdf_open(filenm, /nowrite)
 
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
 
  ; define the date format
  dummy = label_date(date_format=['%D!C%M'])
 
  ; get the number of time elements
  ntime = n_elements(djulian_mod)-1
 
 ; close the NetCDF file
 ncdf_close, nc_file

 ; define plot range
 i_beg = 0
 i_end = ntime-1

 ; make a base plot
 plot, djulian_mod[i_beg:i_end], xrange=[djulian_mod[i_beg],djulian_mod[i_end]], yrange=[0,3], xstyle=9, ystyle=1, $
  xtickformat=['label_date'], xticks=12, ytitle = ytit, xmargin=xmar, ymargin=ymar, xticklen=(-0.02), $
  title=ptitle[ifile]+' ('+cWaterYear+')', /nodata
 plots, [djulian_mod[i_beg],djulian_mod[i_end]], [3,3]

 ; *****
 ; * PLOT THE VALIDATION DATA...
 ; *****************************

 ; open file
 nc_file = ncdf_open(valFile[ifile], /nowrite)

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

  ; get data subset
  caldat, djulian, im, id, iy, ih, imi, asec
  isubset = where( (iy eq iyear[ifile] and im ge 10) or (iy eq iyear[ifile]+1 and im lt 10), nsubset)
  i_beg    = isubset[0]
  i_end    = isubset[nsubset-1]

  ; get the desired variable
  ivar_id = ncdf_varid(nc_file,valName[ifile])
  ncdf_varget, nc_file, ivar_id, valData

  ; remove errors
  if(ifile eq 1)then begin
   if(iyear[ifile] eq 2007) then imissing = where((im ge  7 and im lt 10) or (im eq 6 and id ge 20) or (im eq 10 and id le 5))
   if(iyear[ifile] eq 2008) then imissing = where((im eq 10 and id le 20))
   if(iyear[ifile] le 2008) then valdata[imissing] = -9999.d
  endif

  tvlct, r, g, b, /get
  r[1] = 180
  g[1] = 180
  b[1] = 180
  tvlct, r, g, b

  ; overplot the desired variable
  for itime=i_beg,i_end do begin
   sDepth = valdata[itime]*valMult[ifile]
   if(sDepth gt 0.d)then $
    oplot, [djulian[itime],djulian[itime]], [0.d,sDepth], color=1
  endfor

 ; close the NetCDF file
 ncdf_close, nc_file


 ; *****
 ; * PLOT THE MODEL SIMULATIONS...
 ; *******************************

 ; reset julian day to match model simulations
 djulian = djulian_mod

 ; re-define plot range
 i_beg = 0
 i_end = ntime-1

 ; loop through different parameterizations
 for isuffix=0,n_elements(cSuffix)-1 do begin

  ; declare file
  filenm = file_path+file_pref+cWaterYear+site_name[ifile]+cSuffix[isuffix]+'.nc'

  ; get the snow depth
  nc_file = ncdf_open(filenm, /nowrite)
   ivar_id = ncdf_varid(nc_file,'scalarSnowDepth')
   ncdf_varget, nc_file, ivar_id, scalarSnowDepth
  ncdf_close, nc_file

  ; overplot the snow depth
  oplot, djulian, scalarSnowDepth, color=icolor[iSuffix]

 endfor  ; (looping through experiments)

 ; plot a legend
 if(ifile eq 0) then begin
  xline = [djulian[i_beg]+0.40d*(djulian[i_end]-djulian[i_beg]), $
           djulian[i_beg]+0.55d*(djulian[i_end]-djulian[i_beg])]
  y_top = 2.90
  y_gap = 0.095
  xcolr = [80,250]
  xtext = ['Variable albedo decay','Constant albedo decay']
  for ivar=0,1 do begin
   yline = y_top - float(ivar)*y_gap
   plots, xline, [yline,yline], color=xcolr[ivar]
   xyouts, xline[1]+5., yline-0.025, xtext[ivar], charsize=1.25
  endfor
  ; define the observations
  x1 = xline[0]+10.d
  x2 = xline[1]
  y1 = y_top - float(ivar)*y_gap - y_gap*0.65
  y2 = y_top - float(ivar)*y_gap + y_gap*0.35
  polyfill, [x1,x2,x2,x1], [y1,y1,y2,y2], color=1
  xyouts, x2+5., y1+0.025, 'Observed snow depth', charsize=1.25
 endif

endfor ; loop through sites

; make a figure
write_png, gpath+gname, tvrd(true=1)



stop
end
