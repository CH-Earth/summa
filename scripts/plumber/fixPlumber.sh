#!/opt/local/bin/bash
#
# used to fix some of the plumber files
#
#
# define the new output path
outputPath=/Volumes/d1/mclark/PLUMBER_data/model_output/

# loop through models
for modelName in CHTESSEL COLASSiB.2.0 JULES3.1_altP ORCHIDEE.trunk_r1401; do

 # loop through files
 for outputFile in $( ls  ${outputPath}${modelName}/*.nc ); do

  # compute the net radiation
  ncap2 -A -s 'Rnet=SWnet+LWnet' $outputFile $outputFile
  ncatted -O -a long_name,Rnet,o,c,'Net radiation absorbed within model domain' $outputFile $outputFile 
  echo $outputFile

 done
done
