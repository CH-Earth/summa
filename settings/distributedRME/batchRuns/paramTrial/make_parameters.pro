pro make_parameters

; define the start and end parameter file
iStart = 25
iEnd   = 32

; define the number of files
nFiles = (iEnd - iStart) + 1

; define desired parameters
p01 = create_struct('parName', 'theta_res',          'lowerVal',  0.01,       'upperVal',  0.10,       'defaultVal',   0.139) 
p02 = create_struct('parName', 'theta_sat',          'lowerVal',  0.10,       'upperVal',  0.50,       'defaultVal',   0.550)
p03 = create_struct('parName', 'theta_mp',           'lowerVal',  0.10,       'upperVal',  0.50,       'defaultVal',   0.300)
p04 = create_struct('parName', 'vGn_alpha',          'lowerVal', -1.00,       'upperVal', -0.50,       'defaultVal',  -0.500)
p05 = create_struct('parName', 'vGn_n',              'lowerVal',  1.20,       'upperVal',  1.80,       'defaultVal',   1.300)
p06 = create_struct('parName', 'f_impede',           'lowerVal',  0.00,       'upperVal',  5.00,       'defaultVal',   0.000)
p07 = create_struct('parName', 'k_soil',             'lowerVal',  0.00001,    'upperVal',  0.001,      'defaultVal',   0.0000075)
p08 = create_struct('parName', 'k_macropore',        'lowerVal',  0.001,      'upperVal',  0.100,      'defaultVal',   0.001)
p09 = create_struct('parName', 'kAnisotropic',       'lowerVal',  0.01,       'upperVal', 10.00,       'defaultVal',   1.000)
p10 = create_struct('parName', 'fieldCapacity',      'lowerVal',  0.10,       'upperVal',  0.30,       'defaultVal',   0.200)
p11 = create_struct('parName', 'zScale_TOPMODEL',    'lowerVal',  1.00,       'upperVal',  3.00,       'defaultVal',   2.000)
p12 = create_struct('parName', 'qSurfScale',         'lowerVal', 10.00,       'upperVal',100.00,       'defaultVal',  20.000)
p13 = create_struct('parName', 'critSoilWilting',    'lowerVal',  0.20,       'upperVal',  0.20,       'defaultVal',   0.175)
p14 = create_struct('parName', 'critSoilTranspire',  'lowerVal',  0.25,       'upperVal',  0.25,       'defaultVal',   0.200)
p15 = create_struct('parName', 'heightCanopyTop',    'lowerVal',  0.50,       'upperVal',  0.05,       'defaultVal',   9.500)
p16 = create_struct('parName', 'heightCanopyBottom', 'lowerVal',  0.50,       'upperVal',  0.05,       'defaultVal',   3.000)

; define the number of parameters
nPar = 16

; define the structure
parConstraints = [p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16]
parValues      = fltarr(nPar)

; identify parameter to vary
iSelect = 10

; define the HRU
iHRU=1001

; loop through files
for iFile=iStart,iEnd do begin

 ; define the file name
 fName = 'trial__ex'+strtrim(string(iFile,format='(i2.2)'),2)+'.txt'

 ; open file for writing
 openw, outUnit, fName, /get_lun

 ; write header
 printf, outUnit, '! ***********************************************************************************************************************'
 printf, outUnit, '! ***********************************************************************************************************************'
 printf, outUnit, '! ***** DEFINITION OF TRIAL MODEL PARAMETER VALUES **********************************************************************'
 printf, outUnit, '! ***********************************************************************************************************************'
 printf, outUnit, '! ***********************************************************************************************************************'
 printf, outUnit, '! Note: Lines starting with "!" are treated as comment lines -- there is no limit on the number of comment lines.'
 printf, outUnit, '!'
 printf, outUnit, '! Variable names are important: They must match the variables in the code, and they must occur before the data.'
 printf, outUnit, '!  NOTE: must include information for all HRUs'
 printf, outUnit, '! ***********************************************************************************************************************'

 ; write the parameter names
 printf, outUnit, 'hruIndex', parConstraints[*].ParName, format='(a8,1x,100(a20,1x))'

 ; define the parameter values
 for iPar=0,nPar-1 do begin

  ; populate the parameters
  if(iPar eq iSelect)then begin
   parValues[iPar] = parConstraints[iPar].lowerVal + (parConstraints[iPar].upperVal - parConstraints[iPar].lowerVal)*float(iFile-iStart)/float(nFiles-1)
   print, fName, ': ', parConstraints[iPar].parName, ' = ', parValues[iPar]
  endif else begin
   parValues[iPar] = parConstraints[iPar].defaultVal
  endelse

 endfor  ; looping through parameters

 ; write the parameter values
 printf, outUnit, iHRU, parValues, format='(i8,1x,100(f20.10,1x))'

 ; free up file unit
 free_lun, outUnit


endfor  ; (looping through files)


stop
end
