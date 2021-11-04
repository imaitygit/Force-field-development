#!/bin/bash

#-------------------------------|
# Indrajit Maity                |
# email: i.maity@imperial.ac.uk |
#-------------------------------|

aux_file="/work/e05/e05/imaity/MoTe2_WSe2/mote2_wse2/QE/aux_file"

file1="mote2.in_AA"
file2="wse2.in_AA"
sh="shift_z"

mkdir $sh; cd $sh

# Translate along-z direction
# x = 0; y=0
zmin=7.6; zmax=8.6; dz=0.1
for z in $(seq $zmin $dz $zmax); do
  echo "Two layers are at $z angstrom separation"
  echo "Creating necessary files for QE calculations"

  # create directories
  x=0
  y=0
  mkdir disp_${x}_${y}_${z}
  cd disp_${x}_${y}_${z}
  
  cp $aux_file/*.upf ./
  cp $aux_file/job.slurm ./
  cp $aux_file/$file1 ./
  cp $aux_file/$file2 ./

  # generate the structures
  python $aux_file/auxillary.py "$file1" "$file2" $x $y $z

  # Sumbit jobs
  sbatch job.slurm 
  cd ../

  # Wait for a few minuites before repeating
  sleep 2m
done

