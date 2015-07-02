pro snow_zParamTrial

; used to populate parameter trials for each HRU

; define a line of characters
cLine=''

; define the number of HRUs
nHRU=50

; define desired variables
cDesire = ['hruIndex','winterSAI','summerLAI','theta_res','theta_sat','fieldCapacity','critSoilWilting','critSoilTranspire','vGn_alpha','vGn_n',      'k_soil', 'k_macropore','zScale_TOPMODEL']
xDesire = [         0,       0.0 ,       0.0 ,       0.1 ,     0.375 ,         0.125 ,           0.105 ,             0.300 ,     -0.50 ,   1.5 , 0.0008333333 , 0.0008333333 ,             3.0 ]
nDesire = n_elements(xDesire)

; loop through experiments
for iExp=0,2 do begin

 ; put in the correct parameter sets
 case iExp of

  ; experiment 1
  0: begin
  filenm = 'snow_zParamTrial-exp1.txt'
  xDesire[where(cDesire eq 'k_soil')]          = 0.0000833333
  xDesire[where(cDesire eq 'k_macropore')]     = 0.0000833333
  xDesire[where(cDesire eq 'zScale_TOPMODEL')] = 1.0
  end

  ; experiment 2
  1: begin
  filenm = 'snow_zParamTrial-exp2.txt'
  xDesire[where(cDesire eq 'k_soil')]          = 0.0008333333
  xDesire[where(cDesire eq 'k_macropore')]     = 0.0008333333
  xDesire[where(cDesire eq 'zScale_TOPMODEL')] = 3.0
  end

  ; exfiltration experiment
  2: begin
  filenm = 'snow_zParamTrial-exfiltrate.txt'
  xDesire[where(cDesire eq 'k_soil')]          = 0.00001
  xDesire[where(cDesire eq 'k_macropore')]     = 0.00001
  xDesire[where(cDesire eq 'zScale_TOPMODEL')] = 3.0
  end

  ; check
  else: print, 'cannot find the correct case'

 endcase  ;  put in the correct parameter sets

 ; open file for writing
 openw, outUnit, filenm, /get_lun

 ; write header
 printf, outUnit, '! ***********************************************************************************************************************'
 printf, outUnit, '! ***********************************************************************************************************************'
 printf, outUnit, '! ***** DEFINITION OF TRIAL MODEL PARAMETER VALUES **********************************************************************'
 printf, outUnit, '! ***********************************************************************************************************************'
 printf, outUnit, '! ***********************************************************************************************************************'
 printf, outUnit, '! Note: Lines starting with ! are treated as comment lines -- there is no limit on the number of comment lines.'
 printf, outUnit, '!'
 printf, outUnit, '! Variable names are important: They must match the variables in the code, and they must occur before the data.'
 printf, outUnit, '!  NOTE: must include information for all HRUs'
 printf, outUnit, '! ***********************************************************************************************************************'

 ; write desired variables
 printf, outUnit, cDesire, format='(a8,1x,4(a9,1x),3(a17,1x),2(a9,1x),3(a15,1x))'

 ; loop through HRUs
 for iHRU=1,nHRU do begin
  printf, outUnit, 1000+iHRU, xDesire[1:nDesire-1], format='(i8,1x,4(f9.3,1x),3(f17.3,1x),2(f9.3,1x),2(f15.12,1x),f15.3)'
 endfor


 printf, outUnit, '! ***********************************************************************************************************************'
 printf, outUnit, '! ***********************************************************************************************************************'

 ; free up the file unit
 free_lun, outUnit


endfor  ; looping thru experiments

end
