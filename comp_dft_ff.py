#--------------------------------|
# Author: Indrajit Maity         |
# Email: i.maity@imperial.ac.uk  |
#--------------------------------|

#!/usr/bin/env python3


# To Do
# - Merge x/y/z displacements in a single function 


import sys
import matplotlib.pyplot as plt
import numpy as np
import os


def read_params(p_file):
  """
  Reads the params*.in files and 
  returns the parameters to be tuned.
  IN:
  params_file: file containing parameters
  OUT:
  params: all the parameters, an array
  nv: number of variables
  """
  fp = open(p_file, 'r')
  lines = fp.readlines()
  fp.close()

  for i in range(len(lines)):
    if " variables" in lines[i]:
      nv = eval(lines[i].split()[0])

  try:
    params = np.zeros(nv)
    for i in range(nv):
      params[i] = eval(lines[i+1].split()[0])
  except NameError:
    print("ERROR:")
    print("NameError occured; probably nv is not defined")
    print()
    exit()
  return params, nv


def write_ff_file(type_1, type_2, init, end, params):
  """
  Writes a KC potential file (as required by
  LAMMPS).
  IN:
  type_1: atom type 1, string
  type_1: atom type 2, string
  init: first index in parameters list
  end: the end index in parameters list
  OUT:
  writes out a text file
  """
  # filename based on strings (2-body type)
  ff_file = type_1 + type_2 + str(".KC")
  fp = open(ff_file, 'w')

  # scaling (last term is kept as 1.0)
  fp.write("%s %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"\
           %("#A","B","z0","C0","C2","C4","C","delta",\
            "lambda","A","S"))
  fp.write("%s %s %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %f\n"\
         %(type_1, type_2, params[init],params[init+1],
           params[init+2],params[init+3],params[init+4],
           params[init+5],params[init+6],params[init+7],
           1.0))

  # Following has no effect in calculations
  # It is there to suppress errors in lammps
  if type_1 != type_2:
    fp.write("%s %s %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %f\n"\
            %(type_2, type_1, params[init],params[init+1],
             params[init+2],params[init+3],params[init+4],
             params[init+5],params[init+6],params[init+7],
             0.0))
    fp.write("%s %s %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %f\n"\
            %(type_1, type_1, params[init],params[init+1],
             params[init+2],params[init+3],params[init+4],
             params[init+5],params[init+6],params[init+7],
             0.0))
    fp.write("%s %s %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %f"\
            %(type_2, type_2, params[init],params[init+1],
             params[init+2],params[init+3],params[init+4],
             params[init+5],params[init+6],params[init+7],
             0.0))
  fp.close()
  print()
  print("Writen the file: %s"%(ff_file))
  print("From %d to %d in params.in.*"%(init, end))
  print()


def call_lammps():
  """
  Runs a single-point lammps calculations
  # input file name is: lammps.in 
  """
  fp = open("lammps_command",'r')
  line = fp.readlines()
  fp.close()
  os.system(line[0].strip())
  
 
def gen_lammps_disp_dat(nat, eq_dat, layer2disp, displaced_dat,\
                        dxyz=np.array([0,0,0])):
  """
  IN: 
  nat: number of atoms
  eq_dat: Equilibrium lammps.dat file
          Displacements along x,y,z are carried out
          with respect to eq_dat file
  layer2disp: list/array of atom types to be displaced
  dxyz: displacements along x/y/z directions
            array with default: (0, 0, 0)
  OUT:
  displaced_dat: writes a displaced lammps data file
  """
  # Extract informations from eq_dat file 
  f = open(eq_dat, "r")
  lines = f.readlines()
  f.close()

  # Extract an integer I, till which the eq_dat
  # and output file look identical; see below
  pos = np.zeros((nat, 5), dtype=float)
  for i in range(len(lines)):
    if "Atoms" in lines[i] or "Atoms # atomic" in lines[i]:
      # I is the index
      I = i
      for j in range(i+2, i+2+nat, 1):
        for k in range(5):
          pos[j-i-2][k] = eval(lines[j].split()[k])

  # Output data file
  g = open(displaced_dat, "w")
  
  # write the un-modified lines of lammps data file 
  for i in range(I+2):
    g.write(lines[i])

  # Baseline subtraction
  # Two layers maybe arbitrarily spaced in z-direction
  # Align the interlayer spacing properly
  z1_tmp = []
  z2_tmp = []
  for i in range(nat):
    if int(pos[i][1]) in layer2disp:
      z2_tmp.append(pos[i][4])
    else:
      z1_tmp.append(pos[i][4])
  z1_av = np.average(z1_tmp)
  z2_av = np.average(z2_tmp)
  baseline = z2_av - z1_av

  # write the positions after modifying
  # Note: only atoms belonging to layer2disp
  #       are translated to generate structure
  # BE careful with layer2disp indices!!! 
  # Ensure that indices are counted from 1
  # Otherwise one atom maybe be "badly" displaced.
  for i in range(nat):
    if int(pos[i][1]) in layer2disp:
      if dxyz[2] == 0.0:
        g.write("%d %d %.8f %.8f %.8f\n"%(pos[i][0], pos[i][1],\
              pos[i][2]+dxyz[0], pos[i][3]+dxyz[1], \
              pos[i][4]+dxyz[2]))
      else:
        g.write("%d %d %.8f %.8f %.8f\n"%(pos[i][0], pos[i][1],\
              pos[i][2]+dxyz[0], pos[i][3]+dxyz[1], \
              pos[i][4]+dxyz[2]-(baseline)))
    else:
      g.write("%d %d %.8f %.8f %.8f\n"%(pos[i][0], pos[i][1],\
              pos[i][2], pos[i][3], pos[i][4]))

  g.close()



def gen_lammps_dataset_z(hss,nat,eq_dat,layer2disp,
                         displaced_dat,E_l1,E_l2,\
                         z_min=5.5,z_max=8.6,dz=0.1):
  """
  z_min: minimum of z
  z_max: maximum of z
  dz: steps 
  If not provided defaults are used
  """
  # displacements along z to generate
  z_disp = np.arange(z_min, z_max, dz)

  # compute interlayer binding energy for each configurations
  # with z_disp; hss: high-symmetry-stackings
  if hss == "AA" or hss == "3R":
    print("This is %s high-symmetry stacking"%(hss))
    print("Belongs to 0 degree twist angle")
    print()
  elif hss == "2H" or hss == "BMM" or hss == "BXX":
    print("This is %s high-symmetry stacking"%(hss))
    print("Belongs to 60 degree twist angle")
    print()
  else :
    print("Unrecognized HSS")
    print("Note: this is done intentionally")
    print()
    sys.exit()

  BE = []
  for z in z_disp:

    # displace atoms along z and construct a lammps.dat file
    gen_lammps_disp_dat(nat, eq_dat, layer2disp, displaced_dat,\
                        dxyz=np.array([0,0,z]))

    # run lammps
    call_lammps()
     
    # Extract total energy from lammps output      
    E_bl = extract_toten()

    # E_l1 and E_l2 has to be done on the same
    # single -layer counterparts. It's best to 
    # do it once and be done with it as we only
    # are interested in interlayer energy fitting
    inter_be_patom = (E_bl - E_l1 - E_l2)/nat
    BE.append(inter_be_patom)
  return np.array(BE)


def gen_lammps_dataset_x(hss,nat,eq_dat,layer2disp,
                         displaced_dat,E_l1,E_l2,\
                         x_min=0.0,x_max=4.0,dx=0.1):
  """
  x_min: minimum of x
  x_max: maximum of x
  dx: steps 
  If not provided defaults are used
  """
  # displacements along z to generate
  x_disp = np.arange(x_min, x_max, dx)

  # compute interlayer binding energy for each configurations
  # with x_disp; hss: high-symmetry-stackings
  BE = []
  if hss == "AA" or hss == "3R":
    print("This is %s high-symmetry stacking"%(hss))
    print("Belongs to 0 degree twist angle")
    print()
  elif hss == "2H" or hss == "BMM" or hss == "BXX":
    print("This is %s high-symmetry stacking"%(hss))
    print("Belongs to 60 degree twist angle")
    print()
  else :
    print("Unrecognized HSS")
    print("Note: this is done intentionally")
    print("Neglect this loop if understand\n\
           what you are doing.")
    print()
    sys.exit()

  for x in x_disp:

    # displace atoms along z and construct a lammps.dat file
    gen_lammps_disp_dat(nat, eq_dat, layer2disp, displaced_dat,\
                        dxyz=np.array([x,0,0]))

    # run lammps
    call_lammps()
     
    # Extract total energy from lammps output      
    E_bl = extract_toten()
    
    # E_l1 and E_l2 has to be done on the same
    # single -layer counterparts. It's best to 
    # do it once and be done with it as we only
    # are interested in interlayer energy fitting
    inter_be_patom = (E_bl - E_l1 - E_l2)/nat
    BE.append(inter_be_patom)
  return np.array(BE)


def gen_lammps_dataset_y(hss,nat,eq_dat,layer2disp,
                         displaced_dat,E_l1,E_l2,\
                         y_min=0.0,y_max=4.0,dy=0.1):
  """
  y_min: minimum of y
  y_max: maximum of y
  dy: steps 
  If not provided defaults are used
  """
  # displacements along y to generate
  y_disp = np.arange(y_min, y_max, dy)

  # compute interlayer binding energy for each configurations
  # with y_disp; hss: high-symmetry-stackings
  BE = []
  if hss == "AA" or hss == "3R":
    print("This is %s high-symmetry stacking"%(hss))
    print("Belongs to 0 degree twist angle")
    print()
  elif hss == "2H" or hss == "BMM" or hss == "BXX":
    print("This is %s high-symmetry stacking"%(hss))
    print("Belongs to 60 degree twist angle")
    print()
  else :
    print("Unrecognized HSS")
    print("Note: this is done intentionally")
    print("Neglect this loop if understand\n\
           what you are doing.")
    print()
    sys.exit()

  for y in y_disp:

    # displace atoms along z and construct a lammps.dat file
    gen_lammps_disp_dat(nat, eq_dat, layer2disp, displaced_dat,\
                        dxyz=np.array([0,y,0]))

    # run lammps
    call_lammps()
     
    # Extract total energy from lammps output      
    E_bl = extract_toten()

    # E_l1 and E_l2 has to be done on the same
    # single -layer counterparts. It's best to 
    # do it once and be done with it as we only
    # are interested in interlayer energy fitting
    inter_be_patom = (E_bl - E_l1 - E_l2)/nat
    BE.append(inter_be_patom)
  return np.array(BE)


def extract_toten(lammps_out="log.lammps"):
  """
  IN-
  lammps_out: lammps output file to read
              default: log.lammps which lammps
                       always generates
  OUT-
  toten: total energy at 0-th step
  """

  fp = open(lammps_out, "r")
  lines = fp.readlines()
  fp.close()

  # find/guess the single-point energy
  for i in range(len(lines)):
    # Even if the user customizes output, we expect
    # at-least these two strings:Step and TotEng
    if "Step" in lines[i] and "TotEng" in lines[i]:
      str_id = lines[i].split().index("TotEng")
      return eval(lines[i+1].split()[str_id])


def check_mismatch(arr1, arr2):
  """
  Checks if two lengths are identical or not
  """
  if arr1.size != arr2.size:
    print("ERROR: Check your DFT dataset.")
    sys.exit()

# Calculations
print("----------------------------")
print()
print("    ACTUAL CALCULATIONS     ")
print()
print("----------------------------")


# Individual layer total energy
# MoTe2 single layer from lammps (eV)
E_l1 = -6.23492725699886 

# WSe2 single layer from lammps(eV)
E_l2 = -6.97143986980882 

# Number of atoms
nat = 6  

# Read params file
pfile = sys.argv[1]
params, nv = read_params(pfile)
write_ff_file("Te", "Se", 0, 7, params)
write_ff_file("Mo", "W", 8, 15, params)
write_ff_file("Mo", "Se", 16, 23, params)
write_ff_file("Te", "W", 24, 31, params)

layer2disp = np.array([4, 5, 6])
BE_AA_z = gen_lammps_dataset_z("AA",nat,"lammps.dat_AA",\
                    layer2disp,"lammps.dat_disp",E_l1,E_l2,\
                    z_min=5.5,z_max=8.6,dz=0.1)
BE_3R_z = gen_lammps_dataset_z("3R",nat,"lammps.dat_3R",\
                    layer2disp,"lammps.dat_disp",E_l1,E_l2,\
                    z_min=5.5,z_max=8.6,dz=0.1)
BE_2H_z = gen_lammps_dataset_z("2H",nat,"lammps.dat_2H",\
                    layer2disp,"lammps.dat_disp",E_l1,E_l2,\
                    z_min=5.5,z_max=8.6,dz=0.1)
BE_BMM_z = gen_lammps_dataset_z("BMM",nat,"lammps.dat_BMM",\
                    layer2disp,"lammps.dat_disp",E_l1,E_l2,\
                    z_min=5.5,z_max=8.6,dz=0.1)
BE_BXX_z = gen_lammps_dataset_z("BXX",nat,"lammps.dat_BXX",\
                    layer2disp,"lammps.dat_disp",E_l1,E_l2,\
                    z_min=5.5,z_max=8.6,dz=0.1)
BE_AA_y = gen_lammps_dataset_y("AA",nat,"lammps.dat_AA_z0",\
                    layer2disp,"lammps.dat_disp",E_l1,E_l2,\
                    y_min=0.0,y_max=7.1,dy=0.1)
BE_2H_y = gen_lammps_dataset_y("2H",nat,"lammps.dat_2H_z0",\
                    layer2disp,"lammps.dat_disp",E_l1,E_l2,\
                    y_min=0.0,y_max=7.1,dy=0.1)


def residual(arr1, arr2, weight=1.0, where=[0,31]):
  """
  IN-
  arr1, arr2: DFT and LAMMPS results
  weight: weight to multiply
  where: where to multiply 
  """
  diff = arr1 - arr2
  diff[where[0]:where[1]] = weight*diff[where[0]:where[1]]
  return np.dot(diff, diff)


## Compute the residues
## dft data for identical structures
BE_AA_z_dft = np.loadtxt("DFT/dft_z_AA_be")
check_mismatch(BE_AA_z, BE_AA_z_dft)
BE_3R_z_dft = np.loadtxt("DFT/dft_z_3R_be")
check_mismatch(BE_3R_z, BE_3R_z_dft)
BE_2H_z_dft = np.loadtxt("DFT/dft_z_2H_be")
check_mismatch(BE_2H_z, BE_2H_z_dft)
BE_BXX_z_dft = np.loadtxt("DFT/dft_z_BXX_be")
check_mismatch(BE_BXX_z, BE_BXX_z_dft)
BE_BMM_z_dft = np.loadtxt("DFT/dft_z_BMM_be")
check_mismatch(BE_BMM_z, BE_BMM_z_dft)
BE_AA_y_dft = np.loadtxt("DFT/dft_y_AA_be")
check_mismatch(BE_AA_y, BE_AA_y_dft)
BE_2H_y_dft = np.loadtxt("DFT/dft_y_2H_be")
check_mismatch(BE_2H_y, BE_2H_y_dft)


S_AA_z = residual(BE_AA_z_dft, BE_AA_z, 2, [10, 20])
S_3R_z = residual(BE_3R_z_dft, BE_3R_z, 30, [6, 12])
S_2H_z = residual(BE_2H_z_dft, BE_2H_z, 30, [6, 12])
S_BMM_z = residual(BE_BMM_z_dft, BE_BMM_z, 30, [6, 12])
S_BXX_z = residual(BE_BXX_z_dft, BE_BXX_z, 2, [10, 20])
S_AA_y = residual(BE_AA_y_dft, BE_AA_y, 30, [15, 50])
S_2H_y = residual(BE_2H_y_dft, BE_2H_y, 30, [30, 65])

plt.subplot(2, 2, 1)
plt.plot(BE_AA_z, ls = "--", label="FF(AA)")
plt.plot(BE_3R_z, ls = "--", label=r"FF($\mathrm{B^{Se/Te}}$)")
plt.plot(BE_AA_z_dft, label="DFT(AA)")
plt.plot(BE_3R_z_dft, label=r"DFT($\mathrm{B^{Mo/Se}}$)")
plt.legend(frameon=False)
plt.xlim(5, 20)
plt.ylim(-0.05, 0.05)
#plt.xticks([5, 10, 15, 20], [6, 6.5, 7, 7.5])
plt.subplot(2, 2, 2)
plt.plot(BE_2H_z, ls="--", label="FF(2H)")
plt.plot(BE_BMM_z, ls="--", label=r"FF($\mathrm{B^{Mo/W}}$)")
plt.plot(BE_BXX_z, ls="--", label=r"FF($\mathrm{B^{Se/Te}}$)")
plt.plot(BE_2H_z_dft, label="DFT(2H)")
plt.plot(BE_BMM_z_dft, label=r"DFT($\mathrm{B^{Mo/W}}$)")
plt.plot(BE_BXX_z_dft, label=r"DFT($\mathrm{B^{Se/Te}}$)")
plt.legend(frameon=False)
plt.xlim(5, 20)
plt.ylim(-0.05, 0.05)
#plt.xticks([5, 10, 15, 20], [6, 6.5, 7, 7.5])
plt.subplot(2, 2, 3)
plt.plot(BE_AA_y, ls="--", label="FF(AA)")
plt.plot(BE_AA_y_dft, label="DFT(AA)")
#plt.xticks([0, 10, 20, 30, 40, 50, 60], [0, 1, 2, 3, 4, 5, 6])
plt.legend(frameon=False)
plt.subplot(2, 2, 4)
plt.plot(BE_2H_y, ls="--", label="FF(2H)")
plt.plot(BE_2H_y_dft, label="DFT(2H)")
#plt.xticks([0, 10, 20, 30, 40, 50, 60], [0, 1, 2, 3, 4, 5, 6])
plt.legend(frameon=False)
plt.savefig("fitted.png", dpi=100, bbox_inches="tight", pad_inches=0.1)
#plt.show()

S = S_AA_z + S_3R_z + S_2H_z + S_BMM_z + S_BXX_z + S_AA_y + S_2H_y
#S = S_2H_z + S_AA_z + S_AA_y + S_2H_y
#S = S_2H_y

# Write tot_en in results.out
rfile = sys.argv[2]
fp = open(rfile,'w')
fp.write("%f f"%(S))
fp.close()
