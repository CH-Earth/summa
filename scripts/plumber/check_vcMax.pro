pro check_vcMax

; define plotting parameters
window, 0, xs=1400, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
;!P.MULTI=[0,3,2,0,0]
!P.MULTI=1

; define temperatures
xTemp = dindgen(41)

; define umol2mol conversion factor
uConv = 1000000.d

; define freezing point
Tfreeze = 273.16d

; make a base plot
plot, indgen(5), xrange=[0,40], yrange=[0,100], xstyle=1, ystyle=1, $
 xtitle='Temperature (!eo!nC)', ytitle='vc max', /nodata

; define lookup tables
ixCLM = 1
ixJules = 2

; loop through options
for ixCurve=1,2 do begin

 ; define parameters
 case ixCurve of

  ; define parameters for CLM
  ixCLM: begin
   vcMax25 = [50.d, 60.d, 70., 80.]
   avcmx = 2.4d
   del_S =    710.d
   del_H = 220000.d
   R_gas =      8.314d
   for iPFT=0,3 do begin
    x_num = vcMax25[iPFT]*avcmx^(0.1d*(xTemp - 25.d))
    x_den = 1.d + exp( (-del_H + del_S*(xTemp+Tfreeze)) / (R_gas*(xTemp+Tfreeze)) )
    vcmax = x_num/x_den
    oplot, xTemp, vcMax, color=250
   endfor
  end

  ; define parameters for JULES
  ixJules: begin
   avcmx =  2.0
   eta_e = 0.0008d  ; mol CO2 m-2 s-1 kg C
   eta_0 = [0.046d, 0.033d, 0.073d, 0.060d]   ; [kg C]-1
   vcMax25 = eta_e*eta_0*uConv
   T_low = [ 0., -10.,  0.,  0.]
   T_upp = [36.,  26., 36., 36.]
   for iPFT=0,3 do begin
    x_num = vcMax25[iPFT]*avcmx^(0.1d*(xTemp - 25.d))
    xden1 = (1.d + exp(0.3d*(xTemp - T_upp[iPFT])))
    xden2 = (1.d + exp(0.3d*(T_low[iPFT] - xTemp)))
    vcmax = x_num/(xden1 * xden2)
    oplot, xTemp, vcMax, color=80, linestyle=iPFT, thick=2
   endfor
  end

 endcase

endfor  ; looping through options

; define parameters





end
