pro paramTrial

; *****
; (1) DEFINE PARAMETERS...
; ************************

; define desired parameters
p01 = create_struct('parName', 'theta_res',          'lowerVal',  0.05,       'upperVal',  0.15,       'defaultVal',   0.139)
p02 = create_struct('parName', 'theta_sat',          'lowerVal',  0.45,       'upperVal',  0.60,       'defaultVal',   0.550)
p03 = create_struct('parName', 'theta_mp',           'lowerVal',  0.15,       'upperVal',  0.30,       'defaultVal',   0.300)
p04 = create_struct('parName', 'vGn_alpha',          'lowerVal', -1.00,       'upperVal', -0.50,       'defaultVal',  -0.500)
p05 = create_struct('parName', 'vGn_n',              'lowerVal',  1.20,       'upperVal',  1.80,       'defaultVal',   1.300)
p06 = create_struct('parName', 'f_impede',           'lowerVal',  0.00,       'upperVal',  0.00,       'defaultVal',   0.000)
p07 = create_struct('parName', 'k_soil',             'lowerVal',  0.000001,   'upperVal',  0.00001,    'defaultVal',   0.0000075)
p08 = create_struct('parName', 'k_macropore',        'lowerVal',  0.001,      'upperVal',  0.100,      'defaultVal',   0.001)
p09 = create_struct('parName', 'kAnisotropic',       'lowerVal',  1.00,       'upperVal',  1.00,       'defaultVal',   1.000)
p10 = create_struct('parName', 'fieldCapacity',      'lowerVal',  0.15,       'upperVal',  0.30,       'defaultVal',   0.200)
p11 = create_struct('parName', 'zScale_TOPMODEL',    'lowerVal',  2.00,       'upperVal',  5.00,       'defaultVal',   2.000)
p12 = create_struct('parName', 'qSurfScale',         'lowerVal', 10.00,       'upperVal',100.00,       'defaultVal',  20.000)
p13 = create_struct('parName', 'critSoilWilting',    'lowerVal',  0.20,       'upperVal',  0.20,       'defaultVal',   0.175)
p14 = create_struct('parName', 'critSoilTranspire',  'lowerVal',  0.25,       'upperVal',  0.25,       'defaultVal',   0.200)
p15 = create_struct('parName', 'heightCanopyTop',    'lowerVal',  0.50,       'upperVal',  0.05,       'defaultVal',   9.500)
p16 = create_struct('parName', 'heightCanopyBottom', 'lowerVal',  0.50,       'upperVal',  0.05,       'defaultVal',   3.000)

; get a vector of structures
parConstraints = [p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16]

; identify parameters to vary
cSelect = ['theta_sat','vGn_n','k_soil','k_macropore','zScale_TOPMODEL','qSurfScale']

; define number of parameters to modify
nSelect = n_elements(cSelect)

; define index of parameters to vary
iSelect = intarr(nSelect)

; identify index for each desired parameter
for iDesire=0,nSelect-1 do begin
 ixParam = where(cSelect[iDesire] eq parConstraints[*].parName, nMatch)
 if(nMatch ne 1)then stop, 'unable to identify desired parameter in data structure'
 iSelect[iDesire] = ixParam[0]
endfor

; get default values for selected paramters
xSelect = parConstraints[iSelect].defaultVal

; *****
; (2) DEFINE NEW PARAMETER FILE...
; ********************************

; define a character string
cLine=''

; define the path for the local attributes
attr_path = '/home/mclark/summa/settings/tollgate/localAttributes/'

; define the path for the parameter trials
param_path = '/home/mclark/summa/settings/tollgate/paramTrial/'

; get a list of local attributes files
spawn, 'ls -1 ' + attr_path + 'snow_zLocalAttributes_*.txt', file_list

; get the number of grouped response units
nGRU = n_elements(file_list)

; loop through the basins
for iGRU=0,nGRU-1 do begin

 ; define the name of the local attributes file
 attr_name = file_list[iGRU]

 ; define the name of the parameter trial file
 cFile = strmid(attr_name,strpos(attr_name, '/', /reverse_search)+1)  ; extract the file
 cSuff = strmid(cFile,strpos(cFile, 'basin_', /reverse_search))       ; extract the suffix
 param_name = param_path + 'snow_zParamTrial_' + cSuff
 print, 'param_name = ', param_name

 ; open the local attributes file for reading
 openr, in_unit, attr_name, /get_lun

 ; open the parameter file for writing
 openw, out_unit, param_name, /get_lun

  ; write the parameter header
  printf, out_unit, 'hruIndex', cSelect, format='(a8,1x,100(a20,1x))'

  ; loop through the lines in the veg file
  while ~eof(in_unit) do begin

   ; read a line of data
   readf, in_unit, cLine

   ; split the lines into "words"
   cWords = strsplit(cLine, ' ', /extract)

   ; identify the HRU lines in the attributed file (check that the first word is a number)
   if(stregex(cWords[0], '^[0123456789]',/boolean) eq 1)then begin

    ; get the HRU id (always first)
    ixHRU = long(cWords[0])

    ; write the parameter trial
    printf, out_unit, ixHRU, xSelect, format='(i8,1x,100(f20.10,1x))'

   endif  ; if the word is a number

  endwhile  ; looping through lines in the veg file

 ; free up the file units
 free_lun, in_unit
 free_lun, out_unit

endfor  ; looping through GRUs

stop
end
