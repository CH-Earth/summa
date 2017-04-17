pro compileDecisions

; used to compile a table of decisions from all of the PLUMBER experiments

; define the file path
file_path = '/d1/mclark/PLUMBER_data/model_output/'

; define the file prefix
file_pref = 'SUMMA.1.0.exp.02.'

; define the file suffix
file_suff = 'Fluxnet.1.4.nc'

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

; get number of sites
nSites = n_elements(site_names)

; define desired attributes
model_attr = ['vegeParTbl', $ 
              'soilStress', $ 
              'stomResist', $ 
              'bbTempFunc', $ 
              'bbHumdFunc', $ 
              'bbElecFunc', $ 
              'bbAssimFnc', $ 
              'bbCanIntg8', $ 
              'rootProfil'  ]

; get the number of model attributes
nAtt = n_elements(model_attr)

; define experiments
kExp=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,24,25,26,27,28,29,30,31,32,33]
nExp=n_elements(kExp)

; get an array holding decisions
cDecisions = strarr(nAtt,nExp)

; define a delimiter
cDelim = '--------------------'
cBreak = replicate(cDelim,nAtt)

; loop through sites
for iSite=0,0 do begin

 ; open file for writing
 fileout = 'tableExperiments.txt'
 openw, outUnit, fileout, /get_lun

 ; print decision names
 printf, outUnit, 'ix', model_attr, format='(a2,1x,20(a17,1x))'
 printf, outUnit, '--', cBreak, format='(a2,1x,20(a17,1x))'

 ; loop through experiments
 for iExp=0,nExp-1 do begin

  ; define the experiment
  jExp = kExp[iExp]

  ; define the model name
  model_name = file_pref + strtrim(string(jExp,format='(i3.3)'),2)

  ; define the file name
  file_name = file_path + model_name + '/' + model_name + '_' + site_names[iSite] + file_suff
  ;print, file_name

  ; open the file for reading
  nc_file = ncdf_open(file_name, /nowrite)

   ; loop through model attributes
   for iAttr=0,nAtt-1 do begin

    ; get the model attributes
    ncdf_attget, nc_file, model_attr[iAttr], attVal, /global
    cDecisions[iAttr,iExp] = string(attVal)

   endfor  ; looping through model attributes

  ; close file
  ncdf_close, nc_file

  ; print the table
  printf, outUnit, jExp, cDecisions[*,iExp], format='(i2,1x,20(a17,1x))'

 endfor  ; looping through experiments

 ; free up the file unit
 free_lun, outUnit

endfor  ; looping through sites

spawn, 'cat ' + fileout



stop
end
