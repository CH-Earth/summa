pro check_resFunc

; define plotting parameters
window, 0, xs=1400, ys=800, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=4
!P.COLOR=0
erase, color=255
!P.MULTI=[0,3,2,0,1]

; define named variables for different options
CLM4=1
CLM5=0
Cable=2
Jarvis=3
nModels=4

; define legend
cLegend = strarr(nmodels)
for imodel=0,nmodels-1 do begin
 if(imodel eq CLM4)  then cLegend[imodel] = 'CLM 4'
 if(imodel eq CLM5)  then cLegend[imodel] = 'CLM 5'
 if(imodel eq Cable) then cLegend[imodel] = 'Cable'
 if(imodel eq Jarvis)then cLegend[imodel] = 'Jarvis'
endfor

; define colors
ixColor=[40,90,210,250]
ixLines=[0,4,0,0]

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
joule2umolConv=   4.6d    ; conversion factor from joules to umol photons (umol J-1)
Tfreeze       = 273.16d   ; temerature at freezing (K)
satVPfrz      = 610.8d    ; saturated vapor pressure at freezing (Pa)
wRatio        =   0.622d  ; mixing ratio (-)
Rgas          =   8.314d  ; universal gas constant (J mol-1 K-1)

; define missing value
missingValue = -999.d

; define variables for a given simulation
gb         =        0.05d   ; m s-1
bTran      =        1.0d
sfcprs     =   101325.d

; define default atmospheric o2 and co2 concentration
o2Concentration=0.209d         ; o2 concentration (mol mol-1)  [from Noah-MP  -- also used in Bonan et al., 2011]
;co2Concentration=355.d-6       ; co2 concentration (mol mol-1) [from Noah-MP]
co2Concentration=379.d-6       ; co2 concentration (mol mol-1) [reference value used in Bonan et al., 2011]

; define O2 and CO2 partial pressure
o2  =  o2Concentration * sfcprs  ; 19105.73667200d  ; atmospheric o2 concentration (pa)
co2 = co2Concentration * sfcprs  ;    32.40258630d  ; atmospheric co2 concentration (pa)

; define Michaelis-Menten parameters for CLM4 (Bonan et al., 2011)
Kc25_CLM4  =    30.d    ; Michaelis-Menten constant for CO2 at 25 degrees C (Pa)
Ko25_CLM4  = 30000.d    ; Michaelis-Menten constant for O2 at 25 degrees C (Pa)
qKc_CLM4   =     2.1d   ; factor in the q10 function defining temperature controls on Kc (-)
qKo_CLM4   =     1.2d   ; factor in the q10 function defining temperature controls on Ko (-)

; define Michaelis-Menten parameters for CLM5 (Bonan et al., 2011)
Kc25_CLM5  =   404.90d-6 * sfcprs    ; Michaelis-Menten constant for CO2 at 25 degrees C (Pa)
Ko25_CLM5  =     0.2784d * sfcprs    ; Michaelis-Menten constant for O2 at 25 degrees C (Pa)
Ha_kcCLM5  = 79430.d    ; activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
Ha_koCLM5  = 36380.d    ; activation energy for the Michaelis-Menten constant for O2 (J mol-1)

; define Michaelis-Menten parameters for Cable (Leuning et al. 1995, Table 3)
Kc25_Cable =   302.00d-6 * sfcprs    ; Michaelis-Menten constant for CO2 at 25 degrees C (Pa)
Ko25_Cable =     0.2560d * sfcprs    ; Michaelis-Menten constant for O2 at 25 degrees C (Pa)
Ha_kcCable = 59430.d    ; activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
Ha_koCable = 36000.d    ; activation energy for the Michaelis-Menten constant for O2 (J mol-1)

; define vMax and jmax parameters
vcMax25 =    40.d         ; potential carboxylation rate at 25 degrees C (umol co2 m-2 s-1)
jMax25  =  2.d*vcmax25     ; potential electron transport rate at 25 degrees C (umol co2 m-2 s-1)

; define vcMax parameters from Leuning et al., 1995
; NOTE: these never really worked
;Ha_vcmaxCable = 116300.d    ; activation energy (J mol-1)
;Hd_vcmaxCable = 202900.d    ; deactivation energy (J mol-1)
;Sv_vcmaxCable =    650.d    ; entropy term (J mol-1 K-1)

; define vcMax parameters from Leuning et al., 2002
Ha_vcmaxCable =  73637.d    ; activation energy (J mol-1)
Hd_vcmaxCable = 149252.d    ; deactivation energy (J mol-1)
Sv_vcmaxCable =    486.d    ; entropy term (J mol-1 K-1)

; define jMax parameters from Leuning et al., 2002
Ha_jmaxCable  =  50300.d    ; activation energy (J mol-1)
Hd_jmaxCable  = 152055.d    ; deactivation energy (J mol-1)
Sv_jmaxCable  =    495.d    ; entropy term (J mol-1 K-1)

; define vcmax parameters from Bonan et al., 2011
q10_vcmaxCLM4 =      2.4d   ; factor in the q10 function defining temperature controls on vcmax (-)
Hd_vcmaxCLM4  = 220000.d    ; deactivation energy in high temp inhibition function for vcmax (J mol-1)
Sv_vcmaxCLM4  =    710.d    ; entropy term in high temp inhibition function for vcmax (J K-1 mol-1)

; define vcMax parameters from Bonan et al., 2011
Ha_vcmaxCLM5  =  65330.d    ; activation energy (J mol-1)
Hd_vcmaxCLM5  = 149250.d    ; deactivation energy (J mol-1)
Sv_vcmaxCLM5  =    485.d    ; entropy term (J mol-1 K-1)

; define jMax parameters from Bonan et al. 2011
Ha_jmaxCLM5   =  43540.d    ; activation energy (J mol-1)
Hd_jmaxCLM5   = 152040.d    ; deactivation energy (J mol-1)
Sv_jmaxCLM5   =    495.d    ; entropy term (J mol-1 K-1)

; test
;Ha_vcmaxCable = Ha_vcmaxCLM5  ; activation energy (J mol-1)
;Hd_vcmaxCable = Hd_vcmaxCLM5  ; deactivation energy (J mol-1)
;Sv_vcmaxCable = Sv_vcmaxCLM5  ; entropy term (J mol-1 K-1)

;Kc25_Cable = Kc25_CLM5   ; Michaelis-Menten constant for CO2 at 25 degrees C (Pa)
;Ko25_Cable = Ko25_CLM5   ; Michaelis-Menten constant for O2 at 25 degrees C (Pa)
;Ha_kcCable = Ha_kcCLM5   ; activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
;Ha_koCable = Ha_koCLM5   ; activation energy for the Michaelis-Menten constant for O2 (J mol-1)

; define quantam yield
qYield_CLM4   =      0.05d  ; quantam yield (mol e mol-1 quanta)  [Bonan et al., 2011]
qYield_Cable  =      0.20d  ; quantam yield (mol e mol-1 quanta)  [Table 3, Leuning et al., Plant Cell and Environment 1995]

; define miscellaneous parameters
mp     =      9.d           ; slope parameter
fnf    =      1.d     ; foliage nitrogen factor
gMin   =   2000.d           ; minimum stomatal conductance
e0     =   3500.d           ; scaling factor for humidity function
fJ     =      0.15d         ; parameter to estimate the parameter in the jmax quadratic
Tscale =     10.d           ; scaling factor for q10 temperature function
c_ps2  =      0.7d          ; curvature factor for electron transport (-)

; define smoothing parameters for the assimilation quadratic
theta_cj = 0.98d
theta_ie = 0.95d

; define parameters for the Jarvis scheme
rsmax = 5000.00d ; maximum stomatal resistance (s m-1)
rsmin =  125.00d ; minimum stomatal resistance (s m-1)
topt  =  298.16d ; temperature parameter (K-1)
rgl   =   30.00d ; scaling factor for solar radiation (W-1 m2)
xHS   =   50.00d ; scaling factor for humidity (kg-1 kg)

; define named variables for the type of assimilation
ixRubi = 0
ixLight = 1
ixExport = 2

; define vectors to hold factors and the assimilation
xFac    = dblarr(3)
xAssim  = dblarr(3)

; define named logical variables
yes=0
no=1

; define colimitation
colimitation=yes
;colimitation=no

; define the number of trial values
nTrial = 100

; define the maximum number of iterations
maxiter=10

; define PAR
parVec = 500.d * dindgen(nTrial)/double(ntrial-1)

; define vapor pressure deficit
vpdVec = 3000.d * dindgen(nTrial)/double(ntrial-1)

; define humidity
hsVec = 0.95 * dindgen(nTrial)/double(ntrial-1) + 0.05

; define temperatures
TcVec = 39.d * dindgen(nTrial)/double(ntrial-1) + 0.5d

; define assimilation
asVec = dblarr(nTrial)

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
ixTemp    = 0
ixVPD     = 1
ixPAR     = 2
ixHum     = 3
nx        = 3

; define desired variable for the y-axis
iyAssim   = 0
iyStCond  = 1
iyCarbon  = 2
iyVapPres = 3
ny        = 2

; loop through the x variables
for ixVar=0,nx-1 do begin

 ; loop through the y variables
 for iyVar=0,ny-1 do begin

  ; define number of x tick marks
  if(ixVar eq ixPAR)  then nticks=4 
  if(ixVar eq ixVPD)  then nticks=6
  if(ixVar eq ixHum)  then nticks=4
  if(ixVar eq ixTemp) then nticks=4

  ; define x-axis range
  if(ixVar eq ixPAR)  then xlimits = [0,2000]
  if(ixVar eq ixVPD)  then xlimits = [0,3000]
  if(ixVar eq ixHum)  then xlimits = [0,   1]
  if(ixVar eq ixTemp) then xlimits = [0,  40]

  ; define y-axis range
  if(iyVar eq iyAssim)  then ylimits=[0,  20]
  if(iyVar eq iyStCond) then ylimits=[0,   0.5]
  if(iyVar eq iyCarbon) then ylimits=[0, 100]
  if(iyVar eq iyVapPres)then ylimits=[0,8000]

  ; define x title
  if(iyVar eq ny-1)then begin
   if(ixVar eq ixPAR)  then xVarTitle = 'Absorbed PAR (umol m!e-2!n s!e-1!n)'
   if(ixVar eq ixVPD)  then xVarTitle = 'Vapor pressure deficit (Pa)'
   if(ixVar eq ixHum)  then xVarTitle = 'Canopy air space humidity (fraction)'
   if(ixVar eq ixTemp) then xVarTitle = 'Leaf temperature (!eo!nC)'
  endif else begin
   xVarTitle = ' '
  endelse

  ; define y title
  case ixVar of
   ixTemp: begin
    if(iyVar eq iyAssim)  then yVarTitle='Assimilation!C(mol CO2 m!e-2!n s!e-1!n)'
    if(iyVar eq iyStCond) then yVarTitle='Stomatal condutance!C(mol H2O m!e-2!n s!e-1!n)'
   end
   else: yVarTitle = ' '
  endcase

  ; define x margin
  if(ixVar eq ixTemp)then xMar=[10,-2]
  if(ixVar eq ixHum) then xMar=[ 8, 0]
  if(ixVar eq ixVPD) then xMar=[ 8, 0]
  if(ixVar eq ixPAR) then xMar=[ 6, 2]

  ; define y margin
  if(iyVar eq iyAssim) then yMar=[3,0.5]
  if(iyVar eq iyStCond)then yMar=[4,-0.5]

  ; make a base plot
  plot, indgen(5), xrange=xlimits, yrange=ylimits, xstyle=1, ystyle=1, $
   xtitle=xVarTitle, ytitle=yvarTitle, xticks=nticks, $
   xmargin=xMar, ymargin=yMar, /nodata 

  ; define the axis range
  xr = xlimits[1] - xlimits[0]
  yr = ylimits[1] - ylimits[0]

  ; define the legend
  if(ixVar eq 0)then begin
   x0 = xlimits[0] + 0.05*xr
   x1 = xlimits[0] + 0.25*xr
   for choiceModel=0,nModels-1 do begin
    if(ChoiceModel eq Jarvis and iyVar ne iyStCond)then continue
    y1 = 0.95*yr + ylimits[0] - choiceModel*yr*0.075
    oplot, [x0,x1], [y1,y1], color=ixColor[choiceModel], linestyle=ixLines[choiceModel], thick=3
    xyouts, x1 + 0.01*xr, y1 - 0.025*yr, cLegend[choiceModel], charsize=2
   endfor
  endif

  ; loop through models
  for choiceModel=0,nModels-1 do begin

   ; define Tref
   if(choiceModel eq Cable)then begin
    Tref   =    298.16d ; 25 deg C = reference temperature used in the temperature functions for photosynthesis (K)
   endif else begin
    Tref   =    298.16d ; 25 deg C = reference temperature used in the temperature functions for photosynthesis (K)
   endelse

   ; loop through trial values
   for iTrial=0,nTrial-1 do begin

    ; define humidity
    if(ixVar eq ixHum)then begin
     Tc  =   25.d
     par = 1000.d
     hs  = hsVec[iTrial]
    endif

    ; define vapor pressure deficit
    if(ixVar eq ixVPD)then begin
     Tc  =   25.d
     par = 2000.d / joule2umolConv
     vpd = vpdVec[iTrial]
    endif

    ; define temperature
    if(ixVar eq ixTemp)then begin
     hs  =    1.d
     par = 2000.d / joule2umolConv
     Tc  = TcVec[iTrial]
    endif

    ; define PAR
    if(ixVar eq ixPAR)then begin
     hs  =  1.d
     Tc  = 25.d
     par = parVec[iTrial]
    endif

    ; convert temperatures to kelvin
    Tk = Tc + Tfreeze

    ; define unit conversion (m s-1 --> mol m-2 s-1)
    ; NOTE: Rgas   = J K-1 Mol-1 (J = kg m2 s-2); Tk = K; airpres = Pa (kg m-1 s-2)
    unitConv = sfcprs/(Rgas*Tk)  ; mol m-3

    ; compute the leaf boundary layer resistance (umol-1 m2 s)
    rlb = 1.d/(1000000.d*unitConv*gb)      ; m s-1 --> mol m-2 s-1 --> umol m-2 s-1 --> umol-1 m2 s

    ; compute the saturated vapor pressure (Pa)
    eSat = satVPfrz * exp( (17.27d*Tc)/(237.30d + Tc) )

    ; compute the vapor pressure of the canopy air space (Pa)
    if(ixVar eq ixVPD)then begin
     eAir = eSat - vpd
     hs   = eAir/eSat
     hsVec[iTrial] = hs
    endif else begin
     eAir = hs*eSat
     vpd  = eSat - eAir
    endelse

    ; if CLM4 then constrain humidity
    if(choiceModel eq CLM4 or choiceModel eq CLM5)then hs = max([0.25,hs])

    ; save variables
    eSatVec[iTrial] = eSat
    eAirVec[iTrial] = eAir

    ; compute scaling functions for Jarvis
    if(choiceModel eq Jarvis)then begin

     ; specific humidity
     qAir = wRatio*eAir / (sfcprs - (1.d - wRatio)*eAir)
     qSat = wRatio*eSat / (sfcprs - (1.d - wRatio)*eSat)

     ; mixing ratio
     q2Air = qAir / (1.d + qAir)
     q2Sat = qSat / (1.d + qSat)

     ; contribution due to air temperature
     rct = 1.d - 0.0016d*( (topt - Tk)^2.d )

     ; contribution due to vapor pressure deficit
     rcq = 1.d / (1.d + xHS*(q2Sat - q2Air) )

     ; contribution due to incoming solar radiation
     ff  = 1.d*par/rgl
     rcs = (ff + rsmin/rsmax) / (1.d + ff)

     ; determine canopy resistance (contribution due to all factors)
     rc = rsMin / (rct * rcq * rcs) 

     ; print progress
     print, Tc, vpd, par, eAir, eSat, qAir, qSat, rc, rcs, rct, rcq, unitConv/rc, format=('(2(f9.3,1x),10x,10(f12.7,1x))')

     ; determins stomatal conductance
     gs = 1000000.d*unitConv/rc   ; umol m2 s-1

    endif

    ; compute Michaelis-Menten parameters
    case ChoiceModel of

     CLM4: begin
      Kc = Kc25_CLM4 * call_function('q_T',  Tk, Tref, Tscale, qKc_CLM4)
      Ko = Ko25_CLM4 * call_function('q_T',  Tk, Tref, Tscale, qKo_CLM4)
     end

     CLM5: begin
      Kc = Kc25_CLM5 * call_function('f_T', Tk, Tref, Ha_kcCLM5, Rgas)
      Ko = Ko25_CLM5 * call_function('f_T', Tk, Tref, Ha_koCLM5, Rgas)
     end

     Cable: begin
      Kc = Kc25_Cable * call_function('f_T', Tk, Tref, Ha_kcCable, Rgas)
      Ko = Ko25_Cable * call_function('f_T', Tk, Tref, Ha_koCable, Rgas)
     end

     Jarvis: begin
      Kc = missingValue  ; not used 
      Ko = missingValue  ; not used 
     end

    endcase  ; Michaelis-Menten parameters 

    ; compute vcMax
    case ChoiceModel of

     CLM4: begin
      x_num = call_function('q_T',  Tk, Tref, Tscale, q10_vcmaxCLM4)
      x_den = call_function('fH_T', Tk, Hd_vcmaxCLM4, Sv_vcmaxCLM4, Rgas)
      vcMax = vcMax25*fnf*bTran*x_num/x_den
     end

     CLM5: begin
      f_Tk    = call_function('f_T', Tk, Tref, Ha_vcmaxCLM5, Rgas)
      fH_Tk   = call_function('fH_T', Tk, Hd_vcmaxCLM5, Sv_vcmaxCLM5, Rgas)
      fH_Tref = call_function('fH_T', Tref, Hd_vcmaxCLM5, Sv_vcmaxCLM5, Rgas)
      vcMax   = vcMax25*f_Tk*fH_Tref/fH_Tk
     end

     Cable: begin
      f_Tk    = call_function('f_T', Tk, Tref, Ha_vcmaxCable, Rgas)
      fH_Tk   = call_function('fH_T', Tk, Hd_vcmaxCable, Sv_vcmaxCable, Rgas)
      fH_Tref = call_function('fH_T', Tref, Hd_vcmaxCable, Sv_vcmaxCable, Rgas)
      ;vcMax   = vcMax25*f_Tk/fH_Tk           ; Leuning et al., 1995
      vcMax   = vcMax25*f_Tk*fH_Tref/fH_Tk   ; Leuning et al., 2002
      print, 'vcMax, f_Tk, fH_Tref/fH_Tk, Ha_vcmaxCable, Hd_vcmaxCable, Sv_vcmaxCable = ', $
              vcMax, f_Tk, fH_Tref/fH_Tk, Ha_vcmaxCable, Hd_vcmaxCable, Sv_vcmaxCable, format='(a,1x,10(f12.5,1x))'
     end

     Jarvis: vcMax = -999.  ; not used

    endcase ; vcmax

    ; compute jMax
    case ChoiceModel of

     CLM5: begin
      f_Tk    = call_function('f_T', Tk, Tref, Ha_jmaxCLM5, Rgas)
      fH_Tk   = call_function('fH_T', Tk, Hd_jmaxCLM5, Sv_jmaxCLM5, Rgas)
      fH_Tref = call_function('fH_T', Tref, Hd_jmaxCLM5, Sv_jmaxCLM5, Rgas)
      jMax    = jMax25*f_Tk*fH_Tref/fH_Tk
     end

     Cable: begin
      f_Tk    = call_function('f_T', Tk, Tref, Ha_jmaxCable, Rgas)
      fH_Tk   = call_function('fH_T', Tk, Hd_jmaxCable, Sv_jmaxCable, Rgas)
      fH_Tref = call_function('fH_T', Tref, Hd_jmaxCable, Sv_jmaxCable, Rgas)
      jMax    = jMax25*f_Tk*fH_Tref/fH_Tk
     end

     CLM4:   jMax = -999.  ; not used
     Jarvis: jMax = -999.  ; not used

    endcase  ; jmax

    ; compute electron transport j
    case choiceModel of

     CLM4: begin
      Js = qYield_CLM4*joule2umolConv*PAR  ; available energy: W m-2 -> umol photon m-2 s-1 -> umol e m-2 s-1
     end

     Jarvis: J = missingValue  ; not used

     Cable: begin
      q  = qYield_Cable*joule2umolConv*PAR ; available energy: W m-2 -> umol photon m-2 s-1 -> umol e m-2 s-1
      J  = q*jmax / (q + 2.1d*jmax) 
      Js = J/4.d  ; scaled J
     end

     ; quadratic
     CLM5: begin
      ; The following quadratic equation is used in
      ;   1) Mercade et al., Biogeosciences 2009 (eq A5); and
      ;   2) Bonan et al., JGR 2001 (Table B2).
      ; The equation is also used in
      ;   3) Wittig et al., Global Change Biology 2005 (eq 11)
      ; where they always select the first root (qQuad/aQuad)
      I_ps2 = 0.5d*(1.d - fJ) * joule2umolConv*PAR   ; Farquar (1980), eq 8: PAR absorbed by PS2 (umol photon m-2 s-1)
      ; define coefficients in the quadratic equation
      aQuad = c_ps2            ; quadratic coefficient = cuurvature factor for electron transport
      bQuad = -(I_ps2 + jmax)  ; linear coefficient
      cQuad =  I_ps2 * jmax    ; free term
      ; compute the q term (NOTE: bQuad is always positive)
      bSign = abs(bQuad)/bQuad
      xTemp = bQuad*bQuad - 4.d *aQuad*cQuad
      qQuad = -0.5d * (bQuad + bSign*sqrt(xTemp))
      ; compute roots
      root1 = qQuad / aQuad
      root2 = cQuad / qQuad
      ; select minimum root, required to ensure J=0 when par=0
      ; NOTE: Wittig et al. select the first root, which is the max in all cases I tried
      J  = min([root1,root2])
      Js = J/4.d  ; scaled J
      ; check quadratic
      ;plot_quadratic, aQuad, bQuad, cQuad, root1, root2
      ;print, 'PAR, jmax, root1, root2, max([root1,root2]) = ', $
      ;        PAR, jmax, root1, root2, max([root1,root2]), format='(a,1x,10(f11.4,1x))'
      ;stop
     end  ; CLM5 (quadratic)

    endcase ; options to compute electron transport j

    ; compute the co2 compensation point (pa)
    if(choiceModel eq Cable)then begin
     ;gamma0 = 34.6d-6 * sfcprs     ; mol mol-1 -> Pa
     ;gamma1 =  0.0451d
     ;gamma2 =  0.000347d
     ;cp = gamma0 * (1.d + gamma1*(Tk - Tref) + gamma2*(Tk - Tref)^2.d ) 
     cp = 0.5d*(Kc/Ko)*o2*0.21d
    endif else begin
     cp = 0.5d*(Kc/Ko)*o2*0.21d
    endelse

    ; compute additional control
    awb = Kc*(1.d + o2/ko)
    cp2 = 2.d*cp

    ; define trial value of ci: intercellular co2 (pa)
    if(choiceModel ne Jarvis)then ci = 0.7d*co2

    ; print a marker
    print, '**'
    print, 'vcmax, jmax = ', vcmax, jmax

    ; iterate
    for iter=1,maxiter do begin

     ; don't iterate for Jarvis
     if(choiceModel eq Jarvis)then continue

     ; initialize ci
     ci_old = ci

     ; compute the difference in intercellular co2 compensation from the co2 compensation point
     ciDiff = max([0., ci - cp])

     ; handle constriants in the derivative
     if(ci lt cp)then ciDer=0.d else ciDer=1.d

     ; compute Rubisco-limited assimilation
     xFac[ixRubi]   = vcMax/(ci + awb)
     xAssim[ixRubi] = xFac[ixRubi]*ciDiff

     ; compute light-limited assimilation
     xFac[ixLight]   = Js/(ci + cp2)
     xAssim[ixLight] = xFac[ixLight]*ciDiff

     ; compute export-limited assimilation
     xFac[ixExport]   = 0.5d
     xAssim[ixExport] = xFac[ixExport]*vcmax

     ; no colimitation
     if(colimitation eq no)then begin

      ; compute assimilation
      assim = min(xAssim,ixControl) ; minimum of the three components (ixControl is the index of the minimum)
      x0    = xFac[ixControl]  ; save the factor

      ; compute derivatives in assimilation (no colimitation)
      if(vcmax gt 0.d and Js gt 0.d)then begin
       case ixControl of
        ixRubi:   dA_dc = x0*ciDer - ciDiff*x0*x0/vcmax  ; Rubisco-limited assimilation
        ixLight:  dA_dc = x0*ciDer - ciDiff*x0*x0/Js     ; light-limited assimilation
        ixExport: dA_dc = 0.d                            ; export-limited assimilation
       endcase
      endif else begin
       dA_dc = 0.d
      endelse

     endif else begin

      ; use quadratic function to smooth Rubisco-limited and light-limited (in Sellers et al. 1996; Bonan et al., 2011)
      ; NOTE: originally from Collatz et al. 1991

      ; define ixControl (just in case it is printed)
      ixControl = -999

      ; compute derivatives for individual terms
      if(vcmax gt 0.d)then dAc_dc = xFac[ixRubi]*ciDer - ciDiff*xFac[ixRubi]*xFac[ixRubi]/vcmax else dAc_dc=0.d
      if(Js    gt 0.d)then dAj_dc = xFac[ixLight]*ciDer - ciDiff*xFac[ixLight]*xFac[ixLight]/Js else dAj_dc=0.d
      dAe_dc = 0.d

      ; smooth Rubisco-limitation and light limitation
      quad_smooth, xAssim[ixRubi], xAssim[ixLight], theta_cj, dAc_dc, dAj_dc, xsAssim, dAi_dc, yes
  
      ; smooth intermediate-limitation and export limitation
      quad_smooth, xsAssim, xAssim[ixExport], theta_ie, dAi_dc, dAe_dc, assim, dA_dc, yes
 
      ; check
      ;plot_quadratic, aQuad, bQuad, cQuad, Aroot1, Aroot2
      ;print, xAssim, assim, Aroot1, Aroot2
      ;stop

     endelse

     ; compute co2 concentration at leaf surface (pa)
     x1 = 1.37d * rlb * sfcprs
     cs = co2 - (x1 * assim)

     ; define common term
     if(choiceModel eq Cable)then begin
      css = cs - cp
     endif else begin
      css = cs
     endelse
     x2 = mp*assim*sfcprs/css

     ; define terms for the quadratic function
     if(choiceModel eq CLM4 or choiceModel eq CLM5)then begin
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

     ; compute intercellular co2 concentration (pa)
     x3 = sfcprs*1.65d
     ci = cs - rs*assim*x3

     ; compute derivatives in qquad w.r.t. ci
     dx2_dc = mp*sfcprs*dA_dc*(x1*assim/css + 1.d)/css
     if(choiceModel eq CLM4 or choiceModel eq CLM5)then begin
      dXt_dc = dx2_dc*(rlb*bQuad*2.d - hs*cQuad*4.d)
     endif else begin
      dXt_dc = dx2_dc*(rlb*bQuad*2.d - cQuad*4.d)
     endelse
     dqq_dc = -0.5d * (rlb*dx2_dc + bSign*dXt_dc*0.5d / sqrt(xTemp) )

     ; compute derivatives in rs
     if(root1 gt root2)then begin

      if(choiceModel eq CLM4 or choiceModel eq CLM5)then begin
       drs_dc = (dqq_dc - root1*hs*dx2_dc)/aQuad
       dgs_dc = (hs*dx2_dc - gs*dqq_dc)/qQuad
      endif else begin
       drs_dc = (dqq_dc - root1*dx2_dc)/aQuad
      endelse

     endif else begin

      if(choiceModel eq CLM4 or choiceModel eq CLM5)then begin
       drs_dc = -root2*dqq_dc/qQuad
       dgs_dc = dqq_dc/cQuad
      endif else begin
       drs_dc = -root2*dqq_dc/qQuad
      endelse

     endelse

     ; final derivative
     dci_dc = -dA_dc*(x1 + rs*x3) - x3*assim*drs_dc

     ; check derivatives
     ;func0 = call_function('func_quad', ci_old+0.d, cp, co2, Js, vcMax, theta_cj, awb, rlb, mp, sfcprs, hs, btran, gMin, eair, esat, e0, choiceModel, ixVar)
     ;func1 = call_function('func_quad', ci_old+ dx, cp, co2, Js, vcMax, theta_cj, awb, rlb, mp, sfcprs, hs, btran, gMin, eair, esat, e0, choiceModel, ixVar)
     ;zDeriv = (func1 - func0)/dx
     ;print, 'check derivatives = ', zDeriv, dci_dc
     ;print, 'check derivatives = ', zDeriv, dc_dc

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

     print, 'iter, ixControl, Tc, xAssim, assim, ci_old, co2, ciDiff, xInc = ', $
             iter, ixControl, Tc, xAssim, assim, ci_old, co2, ciDiff, xInc, format='(a,1x,2(i4,1x),20(f14.8,1x))'
 
     ; exit iteration loop
     if(abs(xInc) lt 0.001d)then break

    endfor  ; iterating

    print, 'Tc, assim, gs = ', Tc, assim, gs/1000000.d, format='(a,1x,20(f14.8,1x))'

    ; save assimilation
    asVec[iTrial] = assim

    ; save conductance
    gsVec[iTrial] = gs

    ; compute vapor pressure at the leaf surface (Pa)
    esVec[iTrial] = (eAir*rs + eSat*rlb) / (rlb + rs)

   endfor  ; loop through trial values

   ; only plot stomatal conductance for Jarvis
   if(choiceModel eq Jarvis and iyVar ne iyStCond)then continue

   ; save x var
   if(ixVar eq ixHum)  then xVar=hsVec
   if(ixVar eq ixVPD)  then xVar=vpdVec
   if(ixVar eq ixTemp) then xVar=TcVec
   if(ixVar eq ixPAR)  then xVar=parVec*joule2umolConv

   ; save yvar
   if(iyVar eq iyAssim)   then yVar=asVec
   if(iyVar eq iyStCond)  then yVar=gsVec/1000000.d  ; convert umol --> mol
   if(iyVar eq iyCarbon)  then yVar=ciVec
   if(iyVar eq iyVapPres) then yVar=esVec

   ; plot results
   oplot, xVar, yVar, color=ixColor[choiceModel], linestyle=ixLines[choiceModel], thick=2, min_value=ylimits[0]

  endfor  ; looping through models
  ;if(ixVar eq ixVPD and iyVar gt 0)then stop

  ; plot up carbon at the leaf surface
  if(iyVar eq iyCarbon)then oplot, xVar, csVec, color=160, thick=2, min_value=ylimits[0]

  ; plot up the vapor pressure of the canopy air space and the leaf interior
  if(iyVar eq iyVapPres)then oplot, xVar, eAirVec, color=160, thick=2, min_value=ylimits[0]
  if(iyVar eq iyVapPres)then oplot, xVar, eSatVec, color=210, thick=2, min_value=ylimits[0]

 endfor  ; looping through the y variables

endfor  ; looping through the x variables

write_png, 'figures/stomatalConductance.png', tvrd(true=1)

stop


end

; ****************************************************************************************************************
; ****************************************************************************************************************
; ****************************************************************************************************************
; ****************************************************************************************************************
; ****************************************************************************************************************

function func_quad, ci, cp, co2, J, vcMax, theta_cj, awb, rlb, mp, sfcprs, hs, btran, gmin, eair, esat, e0, choiceModel, ixVar

 ; define named variables for different options
 CLM4=1
 CLM5=0
 Cable=2
 Jarvis=3

 ; define desired variable for the x-axis
 ixTemp    = 0
 ixVPD     = 1
 ixPAR     = 2
 ixHum     = 3

 ; define named variables for the type of assimilation
 ixRubi = 0
 ixLight = 1
 ixExport = 2

 ; define vectors to hold factors and the assimilation
 xFac    = dblarr(3)
 xAssim  = dblarr(3)

 ; define named logical variables
 yes=0
 no=1

 ; define colimitation
 colimitation=yes

 ; define missing value
 valueMissing = -9999.d

 ; compute the difference in intercellular co2 compensation from the co2 compensation point
 ciDiff = max([0., ci - cp])

 ; handle constriants in the derivative
 if(ci lt cp)then ciDer=0.d else ciDer=1.d

 ; compute Rubisco-limited assimilation
 xFac[ixRubi]   = vcMax/(ci + awb)
 xAssim[ixRubi] = xFac[ixRubi]*ciDiff

 ; compute light-limited assimilation
 xFac[ixLight]   = J/(ci + 2.d*cp)
 xAssim[ixLight] = xFac[ixLight]*ciDiff

 ; compute export-limited assimilation
 xFac[ixExport]   = 0.5d
 xAssim[ixExport] = xFac[ixExport]*vcmax

 ; no colimitation
 if(colimitation eq no)then begin

  ; compute assimilation
  assim = min(xAssim,ixControl) ; minimum of the three components (ixControl is the index of the minimum)
  x0    = xFac[ixControl]  ; save the factor

 endif else begin

  ; use quadratic function to smooth Rubisco-limited and light-limited (in Sellers et al. 1996; Bonan et al., 2011)
  ; NOTE: originally from Collatz et al. 1991

  ; smooth Rubisco-limitation and light limitation
  quad_smooth, xAssim[ixRubi], xAssim[ixLight], theta_cj, valueMissing, valueMissing, xsAssim, xJunk, no

  ; smooth intermediate-limitation and export limitation
  quad_smooth, xsAssim, xAssim[ixExport], theta_cj, valueMissing, valueMissing, assim, xJunk, no

 endelse

 ; compute co2 concentration at leaf surface (pa)
 x1 = 1.37d * rlb * sfcprs
 cs = co2 - (x1 * assim)

 ; define common term
 if(choiceModel eq Cable)then begin
  x2 = mp*assim*sfcprs/(cs-cp)
 endif else begin
  x2 = mp*assim*sfcprs/cs
 endelse

 ; terms for the quadratic function
 if(choiceModel eq CLM4 or choiceModel eq CLM5)then begin
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

; **********************************************************************
; **********************************************************************

pro quad_smooth, x1, x2, xsFac, dx1_dc, dx2_dc, xs, dxs_dc, derivDesired

; input
; x1,x2: variables to be smoothed
; xsFac: smoothing factor
; dx1_dc, dx2_dc: derivatives in variables w.r.t. something important
; derivDesired: flag to denote if a derivative is desired

; output
; xs: smoothed variable
; dxs_dc: derivative w.r.t. something important

; procedure
; uses the quadratic of the form
;  xSmooth*xs^2 - (x1 + x2)*xs + x1*x2 = 0

; define named variables
yes=0
no=1

; define the terms in the quadratic
aQuad = xsFac
bQuad = -(x1 + x2)
cQuad = x1*x2

; compute the q term in the quadratic
bSign = abs(bQuad)/bQuad
xTemp = bQuad*bQuad - 4.d *aQuad*cQuad
qquad = -0.5d * (bQuad + bSign*sqrt(xTemp))

; compute roots
root1 = qQuad / aQuad
root2 = cQuad / qQuad
xs    = min([root1,root2])

; check if derivatives are desired
if(derivDesired eq no)then begin

 dxs_dc = 0.d

; compute derivatives
endif else begin

 ; compute derivatives for the terms in the quadratic
 dbq_dc = -(dx1_dc + dx2_dc)
 dcq_dc = x1*dx2_dc + x2*dx1_dc

 ; compute derivatives for xTemp
 dxT_dc = 2.d*(bQuad*dbq_dc) - 4.d*aQuad*dcq_dc
 dqq_dc = -0.5d * (dbq_dc + bsign*dxT_dc/(2.d*sqrt(xTemp)))

 ; compute derivatives in the desired root
 if(root1 lt root2)then begin
  dxs_dc = dqq_dc/aQuad
 endif else begin
  dxs_dc = (dcq_dc - root2*dqq_dc)/qQuad
 endelse

endelse  ; if computing the derivative

end

; **********************************************************************
; **********************************************************************

pro plot_quadratic, aQuad, bQuad, cQuad, root1, root2

; open new graphics window
window, 2, xs=1000, ys=1000, retain=2
device, decomposed=0
!P.MULTI=1

; define scaled x value
n  = 1001
x  = dindgen(n)/double(n)
s1 = -1.e1
s2 =  1.e1
sx = x*(s2 - s1) +s1

; define function
fx = aQuad*sx*sx + bQuad*sx + cQuad

; plot function
plot, sx, fx, /nodata
oplot, sx, fx, color=250
oplot, [s1,s2], [0,0]

; plot roots
plots, root1, 0, psym=sym(6)
plots, root2, 0, psym=sym(6)
xyouts, root1, 0, '1'
xyouts, root2, 0, '2'

end

; **********************************************************************
; **********************************************************************

; q10 function
function q_T, Tk, Tref, Tscale, q10
return, q10^((Tk - Tref)/Tscale)
end

; temperature function
function f_T, Tk, Tref, Ha, Rgas
return, exp( (Ha/(Tref*Rgas)) * (1.d - Tref/Tk) )
end

; high temperature function
function fH_T, Tk, Hd, Sv, Rgas
return, 1.d + exp((Tk*Sv - Hd)/(Tk*Rgas))
end



