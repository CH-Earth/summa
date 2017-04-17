pro check_resFuncSUMMA

; define plotting parameters
window, 0, xs=1400, ys=800, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=4
!P.COLOR=0
erase, color=255
!P.MULTI=[0,3,2,0,1]

; define filepath
file_path = '/home/mclark/summa/output/plumber/stomatalResistance/'

; define columns in the output file
iOutTemp   = 0
iOutVPD    = 1
iOutPAR    = 2
iOutStCond = 3
iOutAssim  = 4
nOutput    = 5

; define named variables for different options
CLM4=1
CLM5=0
Cable=2
nModels=3

; define legend
cLegend = strarr(nmodels)
for imodel=0,nmodels-1 do begin
 if(imodel eq CLM4)  then cLegend[imodel] = 'CLM4'
 if(imodel eq CLM5)  then cLegend[imodel] = 'CLM5'
 if(imodel eq Cable) then cLegend[imodel] = 'Cable'
endfor

; define conversion factor
joule2umolConv=   4.6d    ; conversion factor from joules to umol photons (umol J-1)

; define colors
ixColor=[40,90,210,250]
ixLines=[0,4,0,0]

; define the number of trial values
nTrial = 100

; define desired variable for the x-axis
ixTemp    = 0
ixVPD     = 1
ixPAR     = 2
nx        = 3

; define desired variable for the y-axis
iyAssim   = 0
iyStCond  = 1
ny        = 2

; loop through the x variables
for ixVar=0,nx-1 do begin

 ; loop through the y variables
 for iyVar=0,ny-1 do begin

  ; define number of x tick marks
  if(ixVar eq ixPAR)  then nticks=4 
  if(ixVar eq ixVPD)  then nticks=6
  if(ixVar eq ixTemp) then nticks=4

  ; define x-axis range
  if(ixVar eq ixPAR)  then xlimits = [0,2000]
  if(ixVar eq ixVPD)  then xlimits = [0,3000]
  if(ixVar eq ixTemp) then xlimits = [0,  40]

  ; define y-axis range
  if(iyVar eq iyAssim)  then ylimits=[0,  20]
  if(iyVar eq iyStCond) then ylimits=[0,   0.5]

  ; define x title
  if(iyVar eq ny-1)then begin
   if(ixVar eq ixPAR)  then xVarTitle = 'Absorbed PAR (umol m!e-2!n s!e-1!n)'
   if(ixVar eq ixVPD)  then xVarTitle = 'Vapor pressure deficit (Pa)'
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
    y1 = 0.95*yr + ylimits[0] - choiceModel*yr*0.075
    oplot, [x0,x1], [y1,y1], color=ixColor[choiceModel], linestyle=ixLines[choiceModel], thick=3
    xyouts, x1 + 0.01*xr, y1 - 0.025*yr, cLegend[choiceModel], charsize=2
   endfor
  endif

  ; loop through models
  for choiceModel=0,nModels-1 do begin

   ; define filename
   if(ixVar eq ixTemp)then file_name=cLegend[choiceModel]+'_stomatalResistance.Temp.txt'
   if(ixVar eq ixVPD)then  file_name=cLegend[choiceModel]+'_stomatalResistance.VPD.txt'
   if(ixVar eq ixPAR)then  file_name=cLegend[choiceModel]+'_stomatalResistance.PAR.txt'

   ; define data array
   datarr = fltarr(nOutput,nTrial)

   ; read in model output
   openr, in_unit, file_path+file_name, /get_lun
    readf, in_unit, datarr
   free_lun, in_unit

   ; save x var
   if(ixVar eq ixTemp) then xVar=reform(datarr[iOutTemp,*])
   if(ixVar eq ixVPD)  then xVar=reform(datarr[iOutVPD,*])
   if(ixVar eq ixPAR)  then xVar=reform(datarr[iOutPAR,*])*joule2umolConv

   ; save yvar
   if(iyVar eq iyAssim)   then yVar=reform(datarr[iOutAssim,*])
   if(iyVar eq iyStCond)  then yVar=reform(datarr[iOutStCond,*])

   ; plot results
   oplot, xVar, yVar, color=ixColor[choiceModel], linestyle=ixLines[choiceModel], thick=2, min_value=ylimits[0]

  endfor  ; looping through models

 endfor  ; looping through the y variables

endfor  ; looping through the x variables

write_png, 'figures/stomatalConductanceSUMMA.png', tvrd(true=1)

stop


end
