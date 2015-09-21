pro computeAssimilation

; define plotting parameters
window, 0, xs=1200, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=4
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,3,0,0]

; define legend
cLegend = ['CLM4','Cable']

; define named variables for different options
CLM4 =0
Cable=1

; define colors
ixColor=[80,250]

; define named variables for visualization
doNothing=1000
makeAplot=1001

; define decisions
choiceModel=Cable
visualizeQuadratic=doNothing
;visualizeQuadratic=makeAplot

; define finite difference increment
dx = 1.d-8

; define constants
Tfreeze = 273.16d   ; temerature at freezing (K)
satVPfrz= 610.8d    ; saturated vapor pressure at freezing (Pa)

; define variables for a given simulation
J          =        0.02808771d
hs         =        0.7d
rlb        =        0.00000146d
bTran      =        0.99997730d
vcMax25    =       62.d
sfcprs     =    86143.867d

; define factors to estimate atmospheric o2 and co2 concentration
o2Factor=0.209d         ; empirical factor to obtain partial pressure of o2
co2Factor=355.d-6       ; empirical factor to obtain partial pressure of co2

; define O2 and CO2concentration
o2  =  o2Factor * sfcprs  ; 19105.73667200d  ; atmospheric o2 concentration (pa)
co2 = co2Factor * sfcprs  ;    32.40258630d  ; atmospheric co2 concentration (pa)

; define the vapor pressures
vp_leaf = 2307.72d  ; saturated vapor pressure of the leaf (Pa)
vp_cair =  448.09d  ; saturated vapor pressure of the canopy air space (Pa)
vp_func = (vp_leaf - vp_cair)/1500.d

; define parameters
Kc25  =    30.d     ; Pa
Ko25  = 30000.d     ; Pa
q_Kc  =      2.1d
q_Ko  =      1.2d
avcmx =      2.4d
del_S =    710.d
del_H = 220000.d
R_gas =      8.314d
mp    =      9.d           ; slope parameter
fnf   =      0.666667d     ; foliage nitrogen factor
gMin  =   2000.d           ; minimum stomatal conductance
e0    =   1500.d           ; scaling factor for humidity function

; define the number of trial values
nTrial = 41

; define the maximum number of iterations
maxiter=10

; define humidity
hsVec = 0.9 * dindgen(nTrial)/double(ntrial-1) + 0.1

; define temperatures
TcVec = 30.d * dindgen(nTrial)/double(ntrial-1) + 5.d

; define stomatal conductance
gsVec = dblarr(nTrial)

; define the vapor pressure at the leaf surface
esVec = dblarr(nTrial)

; define the vapor pressure at the leaf interior and the air
eAirVec = dblarr(nTrial)
eSatVec = dblarr(nTrial)

; define the carbon concentrations
ciVec = dblarr(nTrial)
csVec = dblarr(nTrial)

; define desired variable for the x-axis
ixHum     = 1
ixTemp    = 0
nx        = 2

; define desired variable for the y-axis
iyStCond  = 0
iyCarbon  = 1
iyVapPres = 2
ny        = 3

; loop through the x variables
for iyVar=0,ny-1 do begin

 ; loop through the y variables
 for ixVar=0,nx-1 do begin

  ; define x-axis range
  if(ixVar eq ixHum)  then xlimits = [0,1]
  if(ixVar eq ixTemp) then xlimits = [0,40]

  ; define y-axis range
  if(iyVar eq iyStCond) then ylimits=[0,   0.3]
  if(iyVar eq iyCarbon) then ylimits=[0,  50]
  if(iyVar eq iyVapPres)then ylimits=[0,5000]

  ; define x title
  if(ixVar eq ixHum)  then xVarTitle = 'Humidity (fraction)'
  if(ixVar eq ixTemp) then xVarTitle = 'Temperature (!eo!nC)'

  ; define y title
  if(iyVar eq iyStCond) then yvarTitle='Stomatal condutance!C(mol m!e-2!n s!e-1!n)'
  if(iyVar eq iyCarbon) then yvarTitle='Partial pressure of carbon (Pa)'
  if(iyVar eq iyVapPres)then yvarTitle='Vapor pressure (Pa)' 

  ; make a base plot
  plot, indgen(5), xrange=xlimits, yrange=ylimits, xstyle=1, ystyle=1, $
   xtitle=xVarTitle, ytitle=yvarTitle, /nodata 

  ; define the axis range
  xr = xlimits[1] - xlimits[0]
  yr = ylimits[1] - ylimits[0]

  ; define the legend
  x0 = xlimits[0] + 0.05*xr
  x1 = xlimits[0] + 0.25*xr
  for choiceModel=0,1 do begin
   y1 = 0.9*yr + ylimits[0] - choiceModel*yr*0.075
   oplot, [x0,x1], [y1,y1], color=ixColor[choiceModel], thick=3
   xyouts, x1 + 0.01*xr, y1 - 0.025*yr, cLegend[choiceModel], charsize=2
  endfor

  ; loop through models
  for choiceModel=0,1 do begin

   ; loop through trial values
   for iTrial=0,nTrial-1 do begin

    ; define humidity
    if(ixVar eq ixHum)then begin
     Tc = 20.d
     hs = hsVec[iTrial]
    endif

    ; define temperature
    if(ixVar eq ixTemp)then begin
     hs = 0.3d
     Tc = TcVec[iTrial]
    endif

    ; convert temperatures to kelvin
    Tk = Tc + Tfreeze

    ; if CLM4 constrain humidity
    if(choiceModel eq CLM4)then hs = max([0.25,hs])

    ; compute the saturated vapor pressure (Pa)
    eSat = satVPfrz * exp( (17.27d*Tc)/(237.30d + Tc) )

    ; compute the vapor pressure of the canopy air space (Pa)
    eAir = hs*eSat

    ; save variables
    eSatVec[iTrial] = eSat
    eAirVec[iTrial] = eAir

    ; compute Kc and Ko
    Kc = Kc25*q_Kc^(0.1d*(Tc - 25.d))
    Ko = Ko25*q_Ko^(0.1d*(Tc - 25.d))

    ; compute vcMax
    x_num = avcmx^(0.1d*(Tc - 25.d))
    x_den = 1.d + exp( (-del_H + del_S*Tk) / (R_gas*Tk) )
    vcMax = vcMax25*fnf*bTran*x_num/x_den
    ;print, '***'
    ;print, 'Tc, Kc, Ko, vcMax, x_num, x_den = ', Tc, Kc, Ko, vcMax, x_num, x_den

    ; compute additional control
    awb = Kc*(1.d + o2/ko)

    ; compute the co2 compensation point (pa)
    cp = 0.5d*(Kc/Ko)*o2*0.21d

    ; define trial value of ci: intercellular co2 (pa)
    ci = 0.7d*co2

    ; print a marker
    print, '**'

    ; iterate
    for iter=1,maxiter do begin

     ; initialize ci
     ci_old = ci

     ; compute assimilation
     x0 =  vcMax/(ci + awb)
     assim = x0*(ci - cp)
     ;x0 = J/(ci + 2.d * cp)
     ;assim = x0*(ci - cp)
     ;assim = 0.5d * vcmax
     ;print, 'assim = ', assim
 
     ; compute co2 concentration at leaf surface (pa)
     x1 = 1.37d * rlb * sfcprs
     cs = co2 - (x1 * assim)

     ; define common term
     x2 = mp*assim*sfcprs/cs

     ; define terms for the quadratic function
     if(choiceModel eq CLM4)then begin
      aQuad = x2*hs + bTran*gMin
      bQuad = (x2 + bTran*gMin)*rlb - 1.d
      cQuad = -rlb
     endif else begin
      aQuad =  x2 + bTran*gMin*(1.d + (eSat - eAir)/e0)
      bQuad = (x2 + bTran*gMin)*rlb - (eSat - eAir)/e0 - 1.d
      cQuad = -rlb
     endelse
     ;print, 'aQuad, bQuad, cQuad = ', aQuad, bQuad, cQuad, format='(a,1x,10(f20.10,1x))'

     ; compute the q term in the quadratic
     bSign = abs(bQuad)/bQuad
     xTemp = bQuad*bQuad - 4.d *aQuad*cQuad
     qquad = -0.5d * (bQuad + bSign*sqrt(xTemp)) 

     ; compute roots
     root1 = qQuad / aQuad
     root2 = cQuad / qQuad
     rs = max([root1,root2])
 
     ; convert to conductance
     gs = 1.d / rs

     ; save conductance
     gsVec[iTrial] = gs

     ; compute intercellular co2 concentration (pa)
     x2 = sfcprs*1.65d
     ci = cs - rs*assim*x2

     ; compute derivatives in qquad w.r.t. ci
     dA_dc  = x0*x0*(cp + awb)/vcMax
     dx2_dc = mp*sfcprs*dA_dc*(x1*assim/cs + 1.d)/cs
     if(choiceModel eq CLM4)then begin
      dXt_dc = dx2_dc*(rlb*bQuad*2.d - hs*cQuad*4.d)
     endif else begin
      dXt_dc = dx2_dc*(rlb*bQuad*2.d - cQuad*4.d)
     endelse
     dqq_dc = -0.5d * (rlb*dx2_dc + bSign*dXt_dc*0.5d / sqrt(xTemp) )

     ; compute derivatives in rs
     if(root1 gt root2)then begin

      if(choiceModel eq CLM4)then begin
       drs_dc = (dqq_dc - root1*hs*dx2_dc)/aQuad
       dgs_dc = (hs*dx2_dc - gs*dqq_dc)/qQuad
      endif else begin
       drs_dc = (dqq_dc - root1*dx2_dc)/aQuad
      endelse

     endif else begin

      if(choiceModel eq CLM4)then begin
       drs_dc = -root2*dqq_dc/qQuad
       dgs_dc = dqq_dc/cQuad
      endif else begin
       drs_dc = -root2*dqq_dc/qQuad
      endelse

     endelse

     ; final derivative
     dci_dc = -dA_dc*(x1 + rs*x2) - x2*assim*drs_dc

     ; check derivatives
     func0 = call_function('func_quad', ci_old+0.d, cp, co2, J, vcMax, awb, rlb, mp, sfcprs, hs, btran, gMin, eair, esat, e0, choiceModel)
     func1 = call_function('func_quad', ci_old+ dx, cp, co2, J, vcMax, awb, rlb, mp, sfcprs, hs, btran, gMin, eair, esat, e0, choiceModel)
     zDeriv = (func1 - func0)/dx
     print, 'check derivatives = ', zDeriv, dci_dc

     ; visualize roots
     if(visualizeQuadratic eq makeAplot)then begin
      window, 2, xs=1400, ys=1000, retain=2
      device, decomposed=0
      n = 1001
      x = dindgen(n)/double(n)
      s1 = -1.e-4
      s2 =  1.e-4
      sx = x*(s2 - s1) +s1
      fx = aQuad*sx*sx + bQuad*sx + cQuad
      plot, sx, fx, /nodata
      oplot, sx, fx, color=250
      oplot, [s1,s2], [0,0]
      plots, root1, 0, psym=sym(6)
      plots, root2, 0, psym=sym(6)
      xyouts, root1, 0, '1'
      xyouts, root2, 0, '2'
      ;stop
     endif

     ; compute iteration increment
     xInc = (ci - ci_old)/(1.d - dci_dc)

     ; update
     ci = ci_old + xInc

     ; save carbon concentration
     ciVec[iTrial] = ci
     csVec[iTrial] = cs

     print, 'iter, Tc, assim, cs, gsVec[iTrial], ci_old, ci, xInc = ', $
             iter, Tc, assim, cs, gsVec[iTrial], ci_old, ci, xInc, format='(a,1x,i4,1x,20(f16.8,1x))'
 
     ; exit iteration loop
     if(abs(xInc) lt 0.001d)then break

    endfor  ; iterating

    ; compute vapor pressure at the leaf surface (Pa)
    esVec[iTrial] = (eAir*rs + eSat*rlb) / (rlb + rs)

   endfor  ; loop through trial values

   ; save x var
   if(ixVar eq ixHum)  then xVar=hsVec
   if(ixVar eq ixTemp) then xVar=TcVec

   ; save yvar
   if(iyVar eq iyStCond)  then yVar=gsVec/1000000.d  ; convert umol --> mol
   if(iyVar eq iyCarbon)  then yVar=ciVec
   if(iyVar eq iyVapPres) then yVar=esVec

   ; plot results
   oplot, xVar, yVar, color=ixColor[choiceModel], thick=2

  endfor  ; looping through models

  ; plot up carbon at the leaf surface
  if(iyVar eq iyCarbon)then oplot, xVar, csVec, color=160, thick=2

  ; plot up the vapor pressure of the canopy air space and the leaf interior
  if(iyVar eq iyVapPres)then oplot, xVar, eAirVec, color=160, thick=2
  if(iyVar eq iyVapPres)then oplot, xVar, eSatVec, color=210, thick=2

 endfor  ; looping through the y variables

endfor  ; looping through the x variables


stop


end


function func_quad, ci, cp, co2, J, vcMax, awb, rlb, mp, sfcprs, hs, btran, gmin, eair, esat, e0, choiceModel

 ; define named variables for different options
 CLM4 =0
 Cable=1

 ; compute assimilation
 x0 =  vcMax/(ci + awb)
 assim = x0*(ci - cp)
 ;x0 = J/(ci + 2.d * cp)
 ;assim = x0*(ci - cp)
 ;assim = 0.5d * vcmax

 ; compute co2 concentration at leaf surface (pa)
 x1 = 1.37d * rlb * sfcprs
 cs = co2 - (x1 * assim)

 ; define common term
 x2 = mp*assim*sfcprs/cs

 ; terms for the quadratic function
 if(choiceModel eq CLM4)then begin
  aQuad = x2*hs + bTran*gMin
  bQuad = (x2 + bTran*gMin)*rlb - 1.d
  cQuad = -rlb
 endif else begin
  aQuad =  x2 + bTran*gMin*(1.d + (eSat - eAir)/e0)
  bQuad = (x2 + bTran*gMin)*rlb - (eSat - eAir)/e0 - 1.d
  cQuad = -rlb
 endelse

 ; compute q
 bSign = abs(bQuad)/bQuad
 xTemp = bQuad*bQuad - 4.d *aQuad*cQuad
 qquad = -0.5d * (bQuad + bSign*sqrt(xTemp))

 ; compute roots
 root1 = qquad / aQuad
 root2 = cQuad / qQuad

 ; find the maximum root
 rs = max([root1, root2])

 ; compute intercellular co2 concentration (pa)
 x3 = sfcprs*1.65d
 ci = cs - rs*assim*x3

 ; return
 ;return, 1.d / rs
 return, ci

end
