#!/bin/bash
module load StdEnv/2020 intel/2020.1.217 openmpi/4.0.3 cdo/1.9.8

# Define the simulation paths
root_path="/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/concat_summa-"
declare -a method_name=("sundials_1en6" "sundials_1en8") #"be1" "be32" "be64")
arraylength=${#method_name[@]}

# Define a location where we can temporarily store the simulation files so we don't need to swap back and forth
dest_path="/home/avanb/projects/rpp-kshook/avanb/summaWorkflow_data/domain_NorthAmerica/simulation_stats"
mkdir -p $dest_path

# Loop over all simulations and find the mean values we're after
for (( i=0; i<${arraylength}; i++ ));
do
  summa_in="$root_path${method_name[$i]}.nc"
  summa_out="${dest_path}/SUMMA_mean_${method_name[$i]}.nc"

  #echo "index: $i, SUMMA in : $summa_in"
  #echo "index: $i, SUMMA out: $summa_out"

  echo "Working on ${method_name[$i]}"
  cdo timmean -select,name=scalarTotalET,scalarTotalRunoff,gruId,hruId ${summa_in} ${summa_out}
done
