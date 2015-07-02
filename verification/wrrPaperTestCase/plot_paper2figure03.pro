pro plot_paper2figure03

; define plotting parameters
window, 0, xs=1200, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,4,0,0]

; define the color grey
tvlct, r, g, b, /get
r[1] = 180
g[1] = 180
b[1] = 180
tvlct, r, g, b

; define constants
dt      = 3600.d
TFreeze =  273.16d

; define the name of the validation data
valPath = '/home/mclark/test_summa/summa/testCases_data/validationData/'
valFile = valPath + 'ReynoldsCreek_eddyFlux.nc'

; define the path and name of the graphics file
gpath = '/home/mclark/test_summa/summa/verification/zFigures/'
gname = 'Clark_et_al__WRR2015b_figure03.png'

; define file path
file_path = '/home/mclark/test_summa/summa/output/wrrPaperTestCases/'

; define the file prefix
file_pref = ['figure03/vegImpactsWind_2006-2007','figure02/vegImpactsWind_2006-2007']

; define the file suffix
file_suff = ['_riparianAspenExpWindProfile','_riparianAspenWindParamPerturb']

; define HRU
kHRU=[0,2]

; declare file
filenm = file_path + file_pref[0] + file_suff[0] + '.nc'

; *****
; * GET BASIC DATA FROM THE MODEL OUTPUT FILE...
; **********************************************

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
 djulian = bjulian + atime*aoff

 ; get the number of time elements
 ntime = n_elements(djulian)-1

 ; define the date format
 ;dummy = label_date(date_format=['%D-%M!C%H'])
 dummy = label_date(date_format=['%D-%M'])

 i_beg = 0
 i_end = ntime-1

 ; 1 Dec - 31 March
 i_beg = 743
 i_end = 4367

 ; 24 Feb - 21 March
 ;i_beg = 3500
 ;i_end = 4100

 ; define number of desired data points
 ndesire = (i_end - i_beg) + 1

; close NetCDF file
ncdf_close, nc_file

; *****
; * PLOT CUMULATIVE SNOWPACK DRAINAGE...
; **************************************

; plot a "ghost plot" -- just to take care of the first unused space
plot, indgen(5), color=255

; define intensive observation periods
ibeg_iop = [1055,1703]
iend_iop = [1151,1799]

; define colors
icolor=[80,254]

; make a base plot for cumulative melt
plot, djulian[i_beg:i_end], xrange=[djulian[i_beg],djulian[i_end]], yrange=[0,550], xstyle=1, ystyle=1, $
 xtickformat=['label_date'], xticks=5, ytitle = 'Cumulative snowpack drainage (mm)', xmargin=[-56.75,5], ymargin=[-5,1], /nodata

; shade the iops
for iop=0,1 do begin
 x1 = djulian[ibeg_iop[iop]]
 x2 = djulian[iend_iop[iop]]
 y1 = 0
 y2 = 550
 polyfill, [x1,x2,x2,x1], [y1,y1,y2,y2], color=1
 if(iop eq 0)then xyouts, x1 + 0.6d*(x2 - x1), 200., 'First evaluation period', orientation=90, alignment=0.15, charsize=1.5
 if(iop eq 1)then xyouts, x1 + 0.6d*(x2 - x1), 200., 'Second evaluation period', orientation=90, alignment=0.15, charsize=1.5
endfor

; loop through suffixes
for isuff=1,0,-1 do begin

 ; define HRU
 iHRU=kHRU[isuff]

 ; define filename
 filenm = file_path + file_pref[isuff] + file_suff[isuff] + '.nc'

 ; get snowpack drainage
 nc_file = ncdf_open(filenm, /nowrite)
  ivar_id = ncdf_varid(nc_file,'scalarRainPlusMelt')
  ncdf_varget, nc_file, ivar_id, vardata, offset=[iHRU,0], count=[1,ntime]
 ncdf_close, nc_file
 rainPlusMelt = reform(vardata[0,i_beg:i_end])

 ; get rainfall
 nc_file = ncdf_open(filenm, /nowrite)
  ivar_id = ncdf_varid(nc_file,'scalarRainfall')
  ncdf_varget, nc_file, ivar_id, vardata, offset=[iHRU,0], count=[1,ntime]
 ncdf_close, nc_file
 rainfall = reform(vardata[0,i_beg:i_end])

 ; get snow depth
 nc_file = ncdf_open(filenm, /nowrite)
  ivar_id = ncdf_varid(nc_file,'scalarSnowDepth')
  ncdf_varget, nc_file, ivar_id, vardata, offset=[iHRU,0], count=[1,ntime]
 ncdf_close, nc_file
 snowDepth = reform(vardata[0,i_beg:i_end])

 ; get snow melt
 snowMelt = dblarr(n_elements(snowDepth))
 iSnow = where(snowDepth gt 0.001d, nSnow, complement=iBare, nComplement=nBare)
 if(nSnow gt 0)then snowMelt[iSnow] = rainPlusMelt[iSnow]*1000.d
 if(nBare gt 0)then snowMelt[iBare] = rainPlusMelt[iBare]*1000.d - rainfall[iBare]

 ; compute cumulative melt
 cumMelt = dblarr(ndesire)
 cumMelt[0] = snowMelt[0]*dt
 for iCum=1,ndesire-1 do begin
  cumMelt[iCum] = cumMelt[iCum-1] + snowMelt[iCum]*dt
 endfor

 ; plot cumulative melt
 oplot, djulian[i_beg:i_end], cumMelt, color=icolor[isuff], thick=2

endfor

; plot a legend
x0 = djulian[i_beg] + 0.60d*(djulian[i_end] - djulian[i_beg])
x1 = x0 +  1.d
x2 = x1 + 7.d
plots, [x1,x2], [505,505], color=80, thick=2
plots, [x1,x2], [475,475], color=254, thick=2
xyouts, x2+1.d, 470, 'Exponential wind profile extends to the ground surface', charsize=1.25
xyouts, x2+1.d, 500, 'Logarithmic wind profile below vegetation canopy', charsize=1.25


plots, [x0,x0], [450,550]
plots, [x0,djulian[i_end]], [450,450]



; *****
; * PLOT DESIRED VARIABLES...
; ***************************

; define colors
icolor=[80,254]

; define variables
cvar = ['scalarSenHeatGround','scalarSurfaceTemp','scalarSnowDepth']

; define plot titles
ytit = ['Sensible heat (W m!e-2!n)','Surface temperature (!eo!nC)','Snow depth (m)']

; define plot range
ymin = [0,-15,0]
ymax = [1000,5,0.6]

; loop through variables
for ivar=0,n_elements(cvar)-1 do begin

 ; define y-margin
 if(ivar eq 0)then ymar=[-3, 8]
 if(ivar eq 1)then ymar=[ 1, 4]
 if(ivar eq 2)then ymar=[ 5, 0]

 ; loop through intensive observation periods
 for iop=0,1 do begin

  ; define subset
  ibeg_sub = ibeg_iop[iop]
  iend_sub = iend_iop[iop]

  ; make a base plot
  if(ivar eq n_elements(cvar)-1)then begin
   plot, djulian[ibeg_sub:iend_sub], xrange=[djulian[ibeg_sub],djulian[iend_sub]], yrange=[ymin[ivar],ymax[ivar]], xstyle=1, ystyle=1, $
    xtickformat=['label_date'], xticks=4, ytitle = ytit[ivar], xmargin=[10,5], ymargin=ymar, /nodata
  endif else begin
   plot, djulian[ibeg_sub:iend_sub], xrange=[djulian[ibeg_sub],djulian[iend_sub]], yrange=[ymin[ivar],ymax[ivar]], xstyle=1, ystyle=1, $
    xtickname=[' ',' ',' ',' ',' '], xticks=4, ytitle = ytit[ivar], xmargin=[10,5], ymargin=ymar, /nodata
  endelse
  oplot, [djulian[ibeg_sub],djulian[iend_sub]], [0.,0.]

  ; loop through suffixes
  for isuff=0,1 do begin

   ; define HRU
   iHRU=kHRU[isuff]

   ; define filename
   filenm = file_path + file_pref[isuff] + file_suff[isuff] + '.nc'

   ; get data
   nc_file = ncdf_open(filenm, /nowrite)
    ivar_id = ncdf_varid(nc_file,cvar[ivar])
    ncdf_varget, nc_file, ivar_id, vardata, offset=[iHRU,0], count=[1,ntime]
   ncdf_close, nc_file
   vardata = reform(vardata)

   ; scale data
   if(cvar[ivar] eq 'scalarSurfaceTemp')then vardata = vardata - Tfreeze

   ; plot data
   oplot, djulian[i_beg:i_end], vardata[i_beg:i_end], color=icolor[isuff], thick=2

  endfor  ; looping through experiments

  ; plot legend
  if(ivar eq n_elements(cvar)-1)then begin
   xpos = djulian[ibeg_sub] + 0.5d*(djulian[iend_sub] - djulian[ibeg_sub])
   if(iop eq 0)then xyouts, xpos, -0.15, 'First evaluation period', alignment=0.5, charsize=2
   if(iop eq 1)then xyouts, xpos, -0.15, 'Second evaluation period', alignment=0.5, charsize=2
  endif

 endfor  ; looping through iops

endfor  ; looping through variables

; make a figure
write_png, gpath+gname, tvrd(true=1)



stop
end
