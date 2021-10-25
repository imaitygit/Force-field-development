# Before you start running you would require the folllowing:
-----------------------------------------------------------

* Install dakota
* Install LAMMPS
* Prepare DFT dataset
* Prepare high-symmetry stacking file 
  for lammps data run.
* Please go thorugh the comp_dft_data.py
  Some minor things are still done manually.

What the code does?
------------------
It uses DAKOTA to fit the DFT dataset using the parameters.

Running part
-------------
dakota -i kc_moire.in -o kc_moire.out > kc_moire.stdout
