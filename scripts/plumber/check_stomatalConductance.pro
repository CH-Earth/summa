pro check_stomatalConductance

; define plotting parameters
window, 0, xs=700, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=[0,1,2,0,0]

; define file path
fPath = '/d1/mclark/PLUMBER_data/model_output/'

; define the site names
site_names = ['Amplero',     $
              'Blodgett',    $
              'Bugac',       $
              'ElSaler2',    $
              'ElSaler',     $
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

; define the number of sites
nSites  = n_elements(site_names)

; define models
models = ['SUMMA.1.0.exp.02.011', $
          'SUMMA.1.0.exp.02.014' ]

; loop through sites
for iSite=0,nSites-1 do begin
;for iSite=1,1 do begin

 ; set plot window
 window, 1, xs=1400, ys=1200, retain=2
 device, decomposed=0
 LOADCT, 39
 !P.MULTI=[0,2,2,0,0]

 ; read in data
 for ifile=0,1 do begin

  ; define file
  filename = fPath + models[ifile] + '/' + models[ifile] + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'

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

  ; read the intercellular co2 concentration for sunlit leaves
  ;ivar_id = ncdf_varid(nc_file,'scalarIntercellularCO2Sunlit')
  ;ncdf_varget, nc_file, ivar_id, sunlitIntercellularCO2

  ; read the intercellular co2 concentration for shaded leaves
  ;ivar_id = ncdf_varid(nc_file,'scalarIntercellularCO2Shaded')
  ;ncdf_varget, nc_file, ivar_id, shadedIntercellularCO2

  ; read photosynthesis for sunlit leaves
  ivar_id = ncdf_varid(nc_file,'scalarPhotosynthesisSunlit')
  ncdf_varget, nc_file, ivar_id, sunlitPhotosynthesis

  ; read photosynthesis for shaded leaves
  ivar_id = ncdf_varid(nc_file,'scalarPhotosynthesisShaded')
  ncdf_varget, nc_file, ivar_id, shadedPhotosynthesis

  ; read stomatal conductance for sunlit leaves
  ivar_id = ncdf_varid(nc_file,'stomatalConductanceSunlit')
  ncdf_varget, nc_file, ivar_id, sunlitStomatalConductance

  ; read stomatal conductance for shaded leaves
  ivar_id = ncdf_varid(nc_file,'stomatalConductanceShaded')
  ncdf_varget, nc_file, ivar_id, shadedStomatalConductance

  ; read total stomatal conductance
  ivar_id = ncdf_varid(nc_file,'stomatalConductance')
  ncdf_varget, nc_file, ivar_id, totalStomatalConductance

  ; read leaf temperature
  ivar_id = ncdf_varid(nc_file,'scalarCanopyTemp')
  ncdf_varget, nc_file, ivar_id, scalarCanopyTemp

  ; read absorbed solar
  ivar_id = ncdf_varid(nc_file,'scalarCanopyAbsorbedSolar')
  ncdf_varget, nc_file, ivar_id, scalarCanopyAbsorbedSolar

  ; save data
  if(ifile eq 0)then begin
   sunlit_psn01 = reform(sunlitPhotosynthesis)
   shaded_psn01 = reform(shadedPhotosynthesis)
   ;sunlit_c01   = reform(sunlitIntercellularCO2)
   ;shaded_c01   = reform(shadedIntercellularCO2)
   sunlit_g01   = reform(sunlitStomatalConductance)
   shaded_g01   = reform(shadedStomatalConductance)
   total_g01    = reform(totalStomatalConductance)
   swnet_01     = reform(scalarCanopyAbsorbedSolar)
   temp_01      = reform(scalarCanopyTemp)
  endif else begin
   sunlit_psn02 = reform(sunlitPhotosynthesis)
   shaded_psn02 = reform(shadedPhotosynthesis)
   ;sunlit_c02   = reform(sunlitIntercellularCO2)
   ;shaded_c02   = reform(shadedIntercellularCO2)
   sunlit_g02   = reform(sunlitStomatalConductance)
   shaded_g02   = reform(shadedStomatalConductance)
   total_g02    = reform(totalStomatalConductance)
   swnet_02     = reform(scalarCanopyAbsorbedSolar)
   temp_02      = reform(scalarCanopyTemp)
  endelse

 endfor  ; looping through files

 ; define the number of elements to plot
 nplot=ntime
 ;nplot=1000

 ; define yymax
 yymax = 0.01

 ; plot difference in total conductance
 ;plot, total_g01[0:nplot-1], total_g02[0:nplot-1], psym=sym(6), xrange=[0,yymax], yrange=[0,yymax], xstyle=1, ystyle=1, $
 ; xtitle='orig', ytitle='flex'
 ;plots, [0,yymax], [0,yymax], color=80

 ; plot difference in sunlit conductance
 plot, sunlit_g01[0:nplot-1], sunlit_g02[0:nplot-1], psym=sym(6), xrange=[0,yymax], yrange=[0,yymax], xstyle=1, ystyle=1, $
  xtitle='orig', ytitle='flex'
 plots, [0,yymax], [0,yymax], color=80

 ; plot difference in shaded conductance
 plot, shaded_g01[0:nplot-1], shaded_g02[0:nplot-1], psym=sym(6), xrange=[0,yymax], yrange=[0,yymax], xstyle=1, ystyle=1, $
  xtitle='orig', ytitle='flex'
 plots, [0,yymax], [0,yymax], color=80

 ; plot difference in sunlit photosynthesis
 plot, sunlit_g01[0:nplot-1] - sunlit_g02[0:nplot-1], temp_02[0:nplot-1] - 273.16, psym=sym(6), yrange=[0,50], ystyle=1, $ ; xrange=[0,yymax], yrange=[0,yymax], xstyle=1, ystyle=1, $
  xtitle='orig', ytitle='flex'
 plots, [0,0], [0,50], color=80

 ; plot difference in sunlit photosynthesis
 plot, sunlit_g01[0:nplot-1] - sunlit_g02[0:nplot-1], swnet_02[0:nplot-1], psym=sym(6), yrange=[0,1000], ystyle=1, $ ; xrange=[0,yymax], yrange=[0,yymax], xstyle=1, ystyle=1, $
  xtitle='orig', ytitle='flex'
 plots, [0,0], [0,1000], color=80


 ; plot relationship between photosynthesis and conductance
 ;plot, sunlit_psn01[0:nplot-1], sunlit_psn02[0:nplot-1], psym=sym(6), $ ; xrange=[0,yymax], yrange=[0,yymax], xstyle=1, ystyle=1, $
 ; xtitle='orig', ytitle='flex'
 ;plots, [0,0], [0,35], color=80



 ; get dates
 caldat, djulian, im, id, iy, ih, imi, asec

 ; get observation times
 obsTime = ceil(ih*100+imi*1.66666666d)

 ; get the xticks
 xticks = [' ',[strtrim(indgen(7)*3+3,2)],' ']

 ; define ymax
 ymin = 0.00
 ymax = 0.01

 window, 0, xs=1200, ys=1000, retain=2
 device, decomposed=0
 LOADCT, 39
 !P.MULTI=[0,3,4,0,0]

 ; define the months
 cMonth=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

 ; loop through months
 for imonth=1,12 do begin

  ; define the title
  ptitle=site_names[iSite] + ' (' + cMonth[imonth-1] + ')'

  ; make the base plot
  plot, indgen(24)+1, xrange=[0,24], yrange=[ymin,ymax], xstyle=9, ystyle=1, $
    xmargin = [15,2], xticks=8, xtickname=xticks, ytitle='stomatal conductance', $
    xcharsize=1.5, ycharsize=1.5, xticklen=(-0.02), title=ptitle, $
    /nodata
  plots, [0,24], [0,0]
  plots, [0,24], [ymax,ymax]
 
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

 endfor
 stop
endfor  

stop


end
