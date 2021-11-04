#--------------------------------|
# Indrajit Maity                 |
# email: i.maity@imperial.ac.uk  |
#--------------------------------|

import numpy as np
import sys

def read_qe(filename):
  
  f = open(filename, "r")
  lines = f.readlines()
  f.close()

  natom = 3  
  p_c = [[] for j in range(natom)]
  for i in range(len(lines)):
    if "ibrav" in lines[i]:
      a = lines[i].split(',')
      latcon = eval(lines[i+1].split(',')[0].split()[2])
      if eval(a[0].split()[2])== 0:
        print()
        print("The lattice constant for %s is: %.6f"%\
               (filename, latcon))
      else:
        print("Please use ibrav=0 for convenience")
        exit()
    

    elif "ATOMIC_POSITIONS (crystal)" in lines[i]:
      for j in range(natom):
        for k in range(4):
          if k == 0:
            p_c[j].append(str(lines[i+j+1].split()[k]))
          else:
            p_c[j].append(eval(lines[i+j+1].split()[k]))

  return natom, latcon, p_c

def crys2cart(filename):
  """
  IN:
  filename: QE input file
  OUT:
  positions in angstroms
  Easier to compare with LAMMPS this way.
  """
  natom, latcon, p_c = read_qe(filename)

  print()
  print("latcon in Angstrom: %.6f", latcon*0.529177) 
  A = 0.529177 * latcon * np.array([[1.0, 0.0, 0.0],\
                                   [0.5, 0.8660254, 0.0],\
                                   [0.0, 0.0, 7.0362371349]])

  p_ang = np.copy(p_c)
  for i in range(natom):
    for j in range(3):
      p_ang[i][j+1] = (p_c[i][1]*A[0] + p_c[i][2]*A[1] + p_c[i][3]*A[2])[j]

  return natom, latcon, p_ang 

def qe_inp_file(filename1, filename2, dx, dy, z):
  """
  IN-
  filename1: layer1 input file
  filename2: layer2 input file
  dx,dy: displacement along x,y
  z: interlayer spacing between metal layers
  OUT-
  inp_file: input file for QE of bilayer
  """
  natom1, latcon1, p_ang1 = crys2cart(filename1)
  natom2, latcon2, p_ang2 = crys2cart(filename2)

  f = open("pw.in", "w")
  f.write("&CONTROl\n\
  calculation  = 'scf',\n\
  prefix       = 'mote2_wse2',\n\
  pseudo_dir   = './',\n\
  outdir       = './',\n\
/\n\
&SYSTEM\n\
  ibrav     = 0,\n\
  celldm(1) = 6.71422853,\n\
  nat       = %d,\n\
  ntyp      = 4,\n\
  ecutwfc   = 70,\n\
  nspin     = 2\n\
  tot_magnetization = 0.0\n\
  input_dft = 'VDW-DF-C09'\n\
/\n"%(natom1+natom2))

  f.write("&ELECTRONS\n\
  conv_thr    = 1.D-8,\n\
  mixing_beta = 0.7D0,\n\
  electron_maxstep = 500,\n\
  mixing_mode = 'local-TF',\n\
/\n\
&IONS\n\
  ion_dynamics='bfgs',\n\
/\n\
&CELL\n\
 cell_dynamics='bfgs',\n\
 cell_dofree='2Dxy'\n\
 press=0.0,\n\
/\n\
CELL_PARAMETERS {alat}\n\
1.0 0.0 0.0\n\
0.5 0.8660254 0.0\n\
0.0  0.0   7.5301204826\n\
\n\
ATOMIC_SPECIES\n\
  W 183.84   W.upf\n\
  Se 78.96   Se.upf\n\
  Mo 95.95   Mo.upf\n\
  Te 127.6   Te.upf\n\
\n\
ATOMIC_POSITIONS (angstrom)\n")

  # un-shifted
  for i in range(natom1):
    f.write("%s %.6f %.6f %.6f\n"%(str(p_ang1[i][0]),\
              eval(p_ang1[i][1]), eval(p_ang1[i][2]),\
              eval(p_ang1[i][3])))
    z1_av = np.average(p_ang1[:,3].astype(float))
    z2_av = np.average(p_ang2[:,3].astype(float))
    #----------------------------
    disp_z = z - (z2_av - z1_av)
    disp_x = dx; disp_y = dy
    #---------------------------

  for i in range(natom2):
    f.write("%s %.6f %.6f %.6f\n"%(str(p_ang2[i][0]),\
              eval(p_ang2[i][1])+disp_x,\
              eval(p_ang2[i][2])+disp_y,\
              eval(p_ang2[i][3])+disp_z))

  f.write("\nK_POINTS {automatic}\n\
  12 12 1 0 0 0\n")
  f.close()

filename1 = sys.argv[1] 
filename2 = sys.argv[2]
dx = float(sys.argv[3])
dy = float(sys.argv[4])
z  = float(sys.argv[5])

qe_inp_file(filename1, filename2, dx, dy, z)

