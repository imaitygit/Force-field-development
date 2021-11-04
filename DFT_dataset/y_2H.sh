#!/bin/bash

#-------------------------------|
# Indrajit Maity                |
# email: i.maity@imperial.ac.uk |
#-------------------------------|

aux_file="/work/e05/e05/imaity/MoTe2_WSe2/mote2_wse2/QE/aux_file"

file1="mote2.in_2H"
file2="wse2.in_2H"
sh="shift_y"

mkdir $sh; cd $sh

# Translate along-y direction
 x=0
 z=6.85
ymin=4.9; ymax=7.0; dy=0.1
for y in $(seq $ymin $dy $ymax); do
  echo "Two layers are at $z angstrom separation"
  echo "($x $y) in-plane translation"
  echo "Creating necessary files for QE calculations"

  # create directories
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
  sleep 3m
done
