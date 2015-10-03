pro plot_humidityControlSummary

; define plotting parameters
window, 0, xs=1400, ys=900, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=1
erase, color=255
!P.MULTI=[0,7,5,0,0]

; refine colors
tvlct, r, g, b, /get
r[0] = 255
g[0] = 255
b[0] = 255
tvlct, r, g, b

; define file path
fPath = '/d1/mclark/PLUMBER_data/model_output/'

; define the site names
site_names = [$
              'Amplero',     $
              'Blodgett',    $
              ;'Bugac',       $
              ;'ElSaler2',    $
              ;'ElSaler',     $
              ;'Espirra',     $
              ;'FortPeck',    $
              'Harvard',     $
              ;'Hesse',       $
              ;'Howard',      $
              ;'Howlandm',    $
              ;'Hyytiala',    $
              ;'Kruger',      $
              'Loobos',      $
              ;'Merbleue',    $
              ;'Mopane',      $
              ;'Palang',      $
              ;'Sylvania',    $
              ;'Tumba',       $
              ;'UniMich',      $
              ' ']

; define the number of sites
nSites  = n_elements(site_names)-1 ; blank at the end

; define legend text
legendText = ['Linear','Hyperbolic']

; define models
models = ['004', '005' ]

; define the heat map
dx       = 0.15
dy       = 0.75
nx       = 100
ny       = 100
xHeatMap = dx*dindgen(nx) - dx*double(nx)/2.d
yHeatMap = dy*dindgen(ny) - 25.d
zHeatMap = dblarr(nx,ny)

; define x margins
xmar1 = [14,12,10, 8, 6, 15]
xmar2 = [-9,-7,-5,-3,-1,-15]

; define y margins
ymar1 = [ 0, 1, 2, 3, 4, 5, 6]
ymar2 = [ 2, 1, 0,-1,-2,-3,-4]

; loop through sites
for iSite=0,nSites-1 do begin

 ; initialize the number of plots
 kplot=0

 ; define y margins
 ymar=[ymar1[iSite],ymar2[iSite]]

 ; *****
 ; * READ THE DATA...
 ; ******************

 ; read in data
 for ifile=0,1 do begin

  ; define model
  cModel = 'SUMMA.1.0.exp.02.' + models[ifile]

  ; define file
  filename = fPath + cModel + '/' + cModel + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'

  ; open file
  nc_file = ncdf_open(filename, /nowrite)

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

  ; read total stomatal conductance
  ivar_id = ncdf_varid(nc_file,'stomatalConductance')
  ncdf_varget, nc_file, ivar_id, totalStomatalConductance

  ; read leaf temperature
  ivar_id = ncdf_varid(nc_file,'scalarCanopyTemp')
  ncdf_varget, nc_file, ivar_id, scalarCanopyTemp

  ; save data
  if(ifile eq 0)then begin
   total_g01    = reform(totalStomatalConductance)*1000.d  ; m/s -> mm/s
   temp_01      = reform(scalarCanopyTemp)
  endif else begin
   total_g02    = reform(totalStomatalConductance)*1000.d
   temp_02      = reform(scalarCanopyTemp)
  endelse

 endfor  ; looping through files

 ; *****
 ; * COMPUTE THE MEAN DIURNAL CYCLES...
 ; ************************************

 ; get dates
 caldat, djulian, im, id, iy, ih, imi, asec

 ; get observation times
 obsTime = ceil(ih*100+imi*1.66666666d)

 ; get the xticks
 if(iSite eq nSites-1)then xticks = [' ',[strtrim(indgen(7)*3+3,2)],' '] else xticks=replicate(' ', 9)

 ; get the x title
 if(iSite eq nSites-1)then xtitle='Time of day' else xtitle=' '

 ; define ymax
 ymin =  0.00
 ymax = 10.00

 ; define the months
 cMonth=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

 ; loop through months
 for imonth=4,8 do begin

  ; define the x margins
  xmar=[xmar1[kplot],xmar2[kplot]]

  ; define the ytitle
  ytitle=' '
  ;if(kplot eq 0)then ytitle='Stomatal conductance!C(mm/s)' else ytitle=' '

  ; define the number of yticks
  if(kplot eq 0)then yticks=strtrim(string(10.d*dindgen(6)/5.d,format='(f9.0)'),2) else yticks=replicate(' ',6)

  ; make the base plot
  plot, indgen(24)+1, xrange=[0,24], yrange=[ymin,ymax], xstyle=9, ystyle=1, $
    xmargin = xmar, ymargin = ymar, xticks=8, yticks=5,  xtitle=xtitle, ytitle=ytitle, $
    xcharsize=1.5, ycharsize=1.5, xticklen=(-0.02),xtickname=xticks, ytickname=yticks,  $
    /nodata
  plots, [0,24], [0,0]
  plots, [0,24], [ymax,ymax]
  xyouts, 18, 9, cMonth[iMonth-1]
  if(kplot eq 0)then xyouts, 1, 8, site_names[iSite], charsize=3

  ; increment the plot
  kplot=kplot+1

  ; get time of the day
  xtime = dindgen(48)/2.d + 0.25d

  ; get arrays for the mean diurnal cycle
  diurnal_orig = fltarr(49)
  diurnal_new  = fltarr(49)

  ; compute the mean diurnal cycle
  for jhour=0,47 do begin
   ;imatch = where(obsTime eq jhour*50 and im eq imonth and iy eq 2000 and id le 1, nmatch)
   ;imatch = where(obsTime eq jhour*50 and im eq imonth and iy eq 2000, nmatch)
   imatch = where(obsTime eq jhour*50 and im eq imonth, nmatch)
   diurnal_orig[jhour] = total(total_g01[imatch])/float(nmatch)
   diurnal_new[jhour]  = total(total_g02[imatch])/float(nmatch)
  endfor  ; (looping through hours)

  ; wrap-around
  diurnal_orig[48] = diurnal_orig[0]
  diurnal_new[48]  = diurnal_new[0]

  ; plot the diurnal cycle
  oplot, xtime, diurnal_orig,  color=250, thick=3
  oplot, xtime, diurnal_new,  color=80, thick=2, linestyle=4

  ; make a legend
  if(kplot eq 1 and iSite eq 0)then begin
   ytop = 0.7*ymax
   for ifile=0,1 do begin
    ypos = ytop - float(ifile)*ymax*0.1
    if(ifile eq 0)then plots, [1,5], [ypos,ypos], color=250, thick=3
    if(ifile eq 1)then plots, [1,5], [ypos,ypos], color=80, thick=2, linestyle=4
    xyouts, 6, ypos-ymax*0.025, legendText[ifile]
   endfor
  endif

 endfor ; looping through the months

 ; *****
 ; * MAKE A HEAT MAP...
 ; ********************

 ; get difference in stomatal conductance
 xDiff = total_g02 - total_g01

 ; get the mean temperature
 yTemp = (temp_01 + temp_02)/2.d - 273.16d

 ; skip if we have the heat map already
 goto, got_heatmap

 ; get the heat map
 zHeatMap[*,*] = 0.d
 for ix=0,nx-1 do begin
  for iy=0,ny-1 do begin
   iMatch = where(xDiff gt xHeatMap[ix] and xDiff le xHeatMap[ix]+dx and $
                  yTemp gt yHeatMap[iy] and yTemp le yHeatMap[iy]+dy, nMatch)
   zHeatMap[ix,iy] = zHeatMap[ix,iy] + double(nMatch)
  endfor
 endfor

 ; save the heat map
 save, zHeatMap, filename='xIDLsave/humidityControlHeatmap_'+models[0]+'_'+models[1]+'_'+site_names[iSite]+'.sav'

 ; label if have the statistics already
 got_heatmap:

 ; restore the statistics
 restore, 'xIDLsave/humidityControlHeatmap_'+models[0]+'_'+models[1]+'_'+site_names[iSite]+'.sav'

 ; define the x margins
 xmar=[xmar1[kplot],xmar2[kplot]]

 ; define the x title
 if(iSite eq nSites-1)then xtitle='Stomatal conductance!Cdifference(mm/s)' else xtitle=' '

 if(iSite eq nSites-1)then xticks=['-7.5','-2.5','2.5','7.5'] else xticks = replicate(' ', 4)

 ; make a base plot
 plot, xHeatMap, yHeatMap, xrange=[xHeatMap[0],xHeatMap[nx-1]+dx], yrange=[yHeatMap[0],yHeatMap[nx-1]+dy], $
  ;xstyle=1, ystyle=1, xtitle=xtitle, ytitle='Leaf temperature (!eo!nC)', xtickname=xticks, $
  xstyle=1, ystyle=1, xtitle=xtitle, xtickname=xticks, $
  xmargin=xmar, ymargin=ymar, xcharsize=1.5, ycharsize=1.5, xticks=3, /nodata

 ; plot the heat map
 for ix=1,nx-2 do begin
  for iy=1,ny-2 do begin
   xx = [xHeatMap[ix], xHeatMap[ix]+dx, xHeatMap[ix]+dx, xHeatMap[ix]]
   yy = [yHeatMap[iy], yHeatMap[iy], yHeatMap[iy]+dy, yHeatMap[iy]+dy]
   if(zHeatMap[ix,iy] gt 0.)then zz = min([alog10(zHeatMap[ix,iy])*75.d, 250.d]) else zz = 0.d
   polyfill, xx, yy, color=zz
  endfor
 endfor

 ; make a blank plot
 plot, indgen(5), color=!p.background, xmargin=[20,1]


endfor  ; looping through sites

; plot the ytitle
xyouts, 0.025, 0.6, 'Stomatal conductance (mm/s)', orientation=90, alignment=0.5, /normal
xyouts, 0.75, 0.6, 'Leaf temperature (!eo!nC)', orientation=90, alignment=0.5, /normal

; plot a colorbar
ypos1 = 0.24
ypos2 = 0.97
xColors = indgen(251)
xLabels = xColors/75.d
colorbar, 0.93, 0.945, ypos1, ypos2, xLabels, xColors, every=15, charsize=1.5, /nobox, /norm
xyouts, 0.985, 0.5*(ypos1 + ypos2), 'log!i10!n count', alignment=0.5, orientation=90, charsize=4, /normal

; write figure
write_png, 'figures/humidityControlSummary_' + models[0] + '_' + models[1] + '.png', tvrd(true=1)


stop


end
