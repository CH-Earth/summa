pro plot_rootingProfiles

; define plotting parameters
window, 0, xs=2000, ys=1200, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=3
!P.COLOR=0
erase, color=255
!P.MULTI=[0,5,4,0,0]

; define the date format
dummy = label_date(date_format=['%D-%M'])

; define file paths
file_path = '/d1/mclark/PLUMBER_data/model_output/'

; define the model
model_name = 'CABLE_2.0_SLI.vxh599_r553'

; define the site names
site_names = ['Amplero',     $
              'Blodgett',    $
              'Bugac',       $
              'ElSaler',     $
              'ElSaler2',    $
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

site_pfts = ['Grassland',           $       ; 'Amplero',     
             'Evergreen Needleleaf',$       ; 'Blodgett',    
             'Grassland',           $       ; 'Bugac',       
             'Evergreen Needleleaf',$       ; 'ElSaler',    
             'Cropland',            $       ; 'ElSaler2',     
             'Evergreen Broadleaf', $       ; 'Espirra',     
             'Grassland',           $       ; 'FortPeck',    
             'Deciduous Broadleaf', $       ; 'Harvard',     
             'Deciduous Broadleaf', $       ; 'Hesse',       
             'Woody Savanna',       $       ; 'Howard',      
             'Evergreen Needleleaf',$       ; 'Howlandm',    
             'Evergreen Needleleaf',$       ; 'Hyytiala',    
             'Savanna',             $       ; 'Kruger',      
             'Evergreen Needleleaf',$       ; 'Loobos',      
             'Permanent Wetland',   $       ; 'Merbleue',    
             'Woody Savanna',       $       ; 'Mopane',      
             'Evergreen Broadleaf', $       ; 'Palang',      
             'Mixed Forest',        $       ; 'Sylvania',    
             'Evergreen Broadleaf', $       ; 'Tumba',       
             'Deciduous Broadleaf'  ]       ; 'UniMich'      

; define the number of sites
nSites  = n_elements(site_names)

; loop through sites
for iSite=0,nSites-1 do begin

 ; define the file name
 file_name = model_name + '/' + model_name + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'

 ; open files
 ncFileID = ncdf_open(file_path+file_name, /nowrite)

  ; get the depth of each soil layer
  ivar_id = ncdf_varid(ncFileID,'zse')
  ncdf_varget, ncFileID, ivar_id, xVar
  zSoil = reform(xVar)

  ; get the fraction of roots
  ivar_id = ncdf_varid(ncFileID,'froot')
  ncdf_varget, ncFileID, ivar_id, xVar
  fRoot = reform(xVar)

 ; close the netcdf file
 ncdf_close, ncFileID

 ; get the number of soil layers
 nSoil = n_elements(zSoil)

 ; define the title
 plotTitle=site_names[isite] + '!C' + '(' + site_pfts[iSite] + ')'

 ; get arrays for cumulative depth and roots
 cDepth = fltarr(nSoil)
 cRoots = fltarr(nSoil)

 ; get cumulative depth and roots
 cDepth[0] = zSoil[0]
 cRoots[0] = fRoot[0]
 for iLayer=1,nSoil-1 do begin 
  cDepth[iLayer] = cDepth[iLayer-1] + zSoil[iLayer]
  cRoots[iLayer] = cRoots[iLayer-1] + fRoot[iLayer]
 endfor

 ; make a base plot
 plot, cRoots, cDepth, xrange=[0,1], yrange=[2,0], xstyle=1, ystyle=1, $
  xtitle='Fraction of roots', ytitle='Soil Depth', title=plotTitle, $
  ymargin=[4,4], /nodata

 ; plot data
 oplot, cRoots, cDepth, color=250
 oplot, cRoots, cDepth, color=250, psym=sym(6)

 ; define a cumulative rooting profile
 zMax = 1.5d
 zExp = 0.2d
 oplot, (cDepth[0:nSoil-2]/zMax)^zExp, cDepth[0:nSoil-2], color=80
 oplot, (cDepth[0:nSoil-2]/zMax)^zExp, cDepth[0:nSoil-2], color=80, psym=sym(6)

 a = 5.d
 b = 2.d
 d = 2.d
 xx = 1.d - 0.5d*(exp(-a*cDepth[0:nSoil-2]) + exp(-b*cDepth[0:nSoil-2]))
 oplot, xx, cDepth[0:nSoil-2] 

 ; print rooting profiles
 print, '**'
 print, plotTitle
 print, 'zSoil  = ', zSoil,  format='(a,1x,10(f9.3,1x))'
 print, 'fRoot  = ', fRoot,  format='(a,1x,10(f9.3,1x))'
 print, 'cDepth = ', cDepth, format='(a,1x,10(f9.3,1x))'
 print, 'cRoots = ', cRoots, format='(a,1x,10(f9.3,1x))'


endfor  ; looping through sites

; write figure
write_png, 'figures/rootingProfile.png', tvrd(true=1)

stop
end

