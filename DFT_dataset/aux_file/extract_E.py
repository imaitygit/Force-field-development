import numpy as np
import sys

def read_qeout(filename):
  """
  IN-
  filename: QE output file
  OUT-
  E: total energy 
  """
  f = open(filename, "r")
  lines = f.readlines()
  f.close()

  for i in range(len(lines)):
    if "! " in lines[i]:
      try:
        # energy in Ry for bilayer
        E = eval(lines[i].split()[-2])
        return E
      except TypeError:
        print("Problematic file:%s"%(file3))
        exit()

def BE_peratom(file1, file2, file3):
  
  E_l1 = read_qeout(file1)
  E_l2 = read_qeout(file2)
  E = read_qeout(file3)
  # energy in eV
  return (E - (E_l1 + E_l2))*13.60566228/6.

file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]

print(BE_peratom(file1, file2, file3))
