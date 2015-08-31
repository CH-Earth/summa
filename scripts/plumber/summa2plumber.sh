#!/bin/bash
#
# used to convert the summa output files to the plumber format
#

# define named variables to switch between options
ixLoop=1
ixSort=2

# define the desired option to extract the ground heat flux
ixExtract=$ixSort

# define the experiment name
expName=initialPlumberTest

# define the path for the input data
inputPath=/d1/mclark/PLUMBER_data/site_data/met/

# define the output path
outputPath=/home/mclark/summa/output/plumber/

# define the new output path
newOutputPath=/d1/mclark/PLUMBER_data/model_output/

# define the model name
modelName=SUMMA.1.0.exp.01.test

# loop through files in the local attributes folder
for outputFile in $( ls  ${outputPath}orig/*spinup*${expName}.nc ); do

 # split the string using slashes
 IFS='/' read -a strTemp0 <<< "${outputFile}"
 IFS='_' read -a strTemp1 <<< "${strTemp0[-1]}"
 siteName=${strTemp1[0]}
 echo -- $siteName

 # start afresh
 cp ${outputPath}orig/${siteName}*${expName}.nc ${outputPath}temp/

 # loop through desired NetCDF files in the output directory
 for fileName in $( ls  ${outputPath}temp/${siteName}_*${expName}.nc ); do

  # ****************************************************************
  # * extract the ground heat flux
  # ****************************************************************

  # print progress
  echo $fileName

  # 1st option: use looping to extract the ground heat flux
  if [ $ixExtract == $ixLoop ]; then
 
   # remove the Qg variable
   ncks -a -O -x -v Qg $fileName $fileName
 
   # extract the greound heat flux (use a while loop)
   ncap2 -A -v -S extractGroundHeatFlux1.nco $fileName $fileName

  # 2nd option: use sorting to extract the ground heat flux
  else

   # extract the ground heat flux (non-sorted final output)
   ncap2 -A -v -S extractGroundHeatFlux2.nco $fileName $fileName

   # sort the final output
   # the final re-map command does not work within the script
   # (for some reason that I do not understand)
   ncap2 -A -s 'Qg=remap(Qg,srt_map1)' $fileName $fileName

  fi

  # add attributes to the ground heat flux
  ncatted -O -a long_name,Qg,o,c,'Surface ground heat flux' $fileName $fileName
  ncatted -O -a units,Qg,o,c,'W m-2' $fileName $fileName

  # remove undesired variables
  ncks -a -O -x -v nSnow,nLayers,iLayerHeight,iLayerNrgFlux,mLayerTemp,mLayerVolFracLiq,ifcTotoStartIndex,midTotoStartIndex,ifc_qFlux,srt_map1 $fileName $fileName

  # ****************************************************************
  # * compute LAI
  # ****************************************************************

  # compute the total LAI
  ncap2 -A -s 'LAI=scalarLAI+scalarSAI' $fileName $fileName

  # remove undesired variables (scalarLAI and scalarSAI)
  ncks -a -O -x -v scalarLAI,scalarSAI $fileName $fileName

  # ****************************************************************
  # * re-name sensible and latent heat
  # ****************************************************************

  # rename sensible and latent heat
  ncrename -O -v scalarSenHeatTotal,Qh $fileName $fileName
  ncrename -O -v scalarLatHeatTotal,Qle $fileName $fileName

  # ****************************************************************
  # * compute net radiation
  # ****************************************************************

  # compute the net shortwave and the net longwave radiation
  ncap2 -A -s 'SWnet=scalarCanopyAbsorbedSolar+scalarGroundAbsorbedSolar' $fileName $fileName
  ncap2 -A -s 'LWnet=scalarLWNetCanopy+scalarLWNetGround' $fileName $fileName

  # compute the net radiation
  ncap2 -A -s 'Rnet=scalarCanopyAbsorbedSolar+scalarGroundAbsorbedSolar+scalarLWNetCanopy+scalarLWNetGround' $fileName $fileName

  # modify the long names
  ncatted -O -a long_name,Rnet,o,c,'Net radiation absorbed within model domain' $fileName $fileName
  ncatted -O -a long_name,SWnet,o,c,'SW radiation absorbed within model domain' $fileName $fileName
  ncatted -O -a long_name,LWnet,o,c,'LW radiation absorbed within model domain' $fileName $fileName

 done

 # ****************************************************************
 # * concatanate the files
 # ****************************************************************

 # concatanate
 ncrcat -O ${outputPath}temp/${siteName}_*spinup*.nc ${outputPath}temp/${siteName}_*-*.nc ${outputPath}SUMMA.1.0_${siteName}Fluxnet.1.4.nc

done
