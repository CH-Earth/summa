pro plot_humidityControlScatter

; define plotting parameters
window, 1, xs=2000, ys=1200, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=1
erase, color=255
!P.MULTI=[0,5,4,0,0]

; refine colors
tvlct, r, g, b, /get
r[0] = 255
g[0] = 255
b[0] = 255
tvlct, r, g, b

; define file path
fPath = '/d1/mclark/PLUMBER_data/model_output/'

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

; define the number of sites
nSites  = n_elements(site_names)

; define models
models = ['002', $
          '003' ]

; define the heat map
dx       = 0.15
dy       = 0.75
nx       = 100
ny       = 100
xHeatMap = dx*dindgen(nx) - dx*double(nx)/2.d
yHeatMap = dy*dindgen(ny) - 25.d
zHeatMap = dblarr(nx,ny)

; loop through sites
for iSite=0,nSites-1 do begin

 ; skip if we have the heat map already
 ;goto, got_heatmap

 ; initialize the heat map
 zHeatMap[*,*] = 0.d

 ; read in data
 for ifile=0,1 do begin

  ; define model
  cModel = 'SUMMA.1.0.exp.02.' + models[ifile]

  ; define file
  filename = fPath + cModel + '/' + cModel + '_' + site_names[iSite] + 'Fluxnet.1.4.nc'

  ; open file
  nc_file = ncdf_open(filename, /nowrite)

  ; read total stomatal conductance
  ivar_id = ncdf_varid(nc_file,'stomatalConductance')
  ncdf_varget, nc_file, ivar_id, totalStomatalConductance

  ; read leaf temperature
  ivar_id = ncdf_varid(nc_file,'scalarCanopyTemp')
  ncdf_varget, nc_file, ivar_id, scalarCanopyTemp

  ; save data
  if(ifile eq 0)then begin
   total_g01    = reform(totalStomatalConductance)*1000.d
   temp_01      = reform(scalarCanopyTemp)
  endif else begin
   total_g02    = reform(totalStomatalConductance)*1000.d
   temp_02      = reform(scalarCanopyTemp)
  endelse

 endfor  ; looping through files

 ; get difference in stomatal conductance
 xDiff = total_g02 - total_g01

 ; get the mean temperature
 yTemp = (temp_01 + temp_02)/2.d - 273.16d

 ; get the heat map
 for ix=0,nx-1 do begin
  for iy=0,ny-1 do begin
   iMatch = where(xDiff gt xHeatMap[ix] and xDiff le xHeatMap[ix]+dx and $
                  yTemp gt yHeatMap[iy] and yTemp le yHeatMap[iy]+dy, nMatch)
   zHeatMap[ix,iy] = zHeatMap[ix,iy] + double(nMatch)
  endfor
 endfor

 ; save the heat map
 save, zHeatMap, filename='xIDLsave/humidityControlHeatmap_'+models[0]+'_'+models[1]+'_'+site_names[iSite]+'.sav'

 ; label if have the statistics already
 got_heatmap:

 ; restore the statistics
 restore, 'xIDLsave/humidityControlHeatmap_'+models[0]+'_'+models[1]+'_'+site_names[iSite]+'.sav'

 ; make a base plot
 plot, xHeatMap, yHeatMap, xrange=[xHeatMap[0],xHeatMap[nx-1]+dx], yrange=[yHeatMap[0],yHeatMap[nx-1]+dy], $
  xstyle=1, ystyle=1, xtitle='Difference in stomatal conductance', ytitle='Temperature', $
  title=site_names[isite], /nodata

 ; plot the heat map
 for ix=1,nx-2 do begin
  for iy=1,ny-2 do begin
   xx = [xHeatMap[ix], xHeatMap[ix]+dx, xHeatMap[ix]+dx, xHeatMap[ix]]
   yy = [yHeatMap[iy], yHeatMap[iy], yHeatMap[iy]+dy, yHeatMap[iy]+dy]
   if(zHeatMap[ix,iy] gt 0.)then zz = min([alog10(zHeatMap[ix,iy])*75.d, 250.d]) else zz = 0.d
   polyfill, xx, yy, color=zz
  endfor
 endfor

endfor  ; looping through sites  

; write figure
write_png, 'figures/humidityControlHeatMap_' + models[0] + '_' + models[1] + '.png', tvrd(true=1)


stop


end
