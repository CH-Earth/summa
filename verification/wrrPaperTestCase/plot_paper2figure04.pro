pro plot_paper2figure04

; define plotting parameters
window, 1, xs=800, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=[0,1,2]

; define temperatures
tair = -dindgen(61)/5.d + 2.d

; define density
dens = 67.92d + 51.25d*exp(Tair/2.59d)

; define temperature impact on interception capacity
delI = 0.27d + 46.d/dens

; plot temperature - temp func
plot, tair, delI, xrange=[-10,2], yrange=[0,1], xstyle=1, ystyle=9, $
 xtitle='Temperature (!eo!nC)', ytitle='Scaled interception capacity (-)', $
 xmargin=[10,10], title='Hedstrom and Pomeroy (1988)', /nodata

oplot, tair, delI, color=250, thick=2

; make a second axis for snow density
axis, yaxis=1, ytitle='Snow density (kg m!e-3!n)', yrange=[0,200], ystyle=1, color=80, /save
oplot, tair, dens, color=80, thick=2

; define leaf area ratio
lrat = dblarr(61)
lrat[where(tair gt -1.d)] = 4.d
lrat[where(tair le -3.d)] = 1.0d
ipos = where(tair gt -3.d and tair le -1.d)
lrat[ipos] = 1.5d*tair[ipos] + 5.5d

; plot temperature - VIC temp func
plot, tair, lrat, xrange=[-10,2], yrange=[1,4], xstyle=1, ystyle=1, $
 xtitle='Temperature (!eo!nC)', ytitle='Scaled interception capacity (-)', $
 xmargin=[10,10], title='Andreadis et al. (2009)', /nodata

oplot, tair, lrat, color=250, thick=2

; make a figure
write_png, 'zFigures/Clark_et_al__WRR2015b_figure04.png', tvrd(true=1)

stop
end
