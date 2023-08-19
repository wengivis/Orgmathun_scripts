#!/usr/bin/env python
#
# Obtain the HOMA and HOMER value for the user defined ring in xyz geometry
#

#
# Import libraries
#

import math
import sys
import os

#
# Definition of HOMA and HOMER parameters
#

# Definition of HOMA parameters from Chem. Rev., 2014, 114, 6383-6422

opt_cc_homa = 1.388
a_cc_homa = 257.700

opt_cn_homa = 1.334
a_cn_homa = 93.520

opt_co_homa = 1.265
a_co_homa = 157.380

opt_cb_homa = 1.4235
a_cb_homa = 104.507

opt_cs_homa = 1.677
a_cs_homa = 94.090

opt_cp_homa = 1.698
a_cp_homa = 118.910

opt_nn_homa = 1.309
a_nn_homa = 130.330

opt_nb_homa = 1.402
a_nb_homa = 72.030

opt_no_homa = 1.248
a_no_homa = 57.210

# Definition of HOMER parameters from Phys. Chem. Chem. Phys., 2023, 25, 16763-16771

opt_cc_homer = 1.437
a_cc_homer = 950.74

opt_cn_homer = 1.390
a_cn_homer = 506.43

opt_nn_homer = 1.375
a_nn_homer = 187.36

opt_co_homer = 1.379
a_co_homer = 164.96


#
# Definition of functions.
#

# Function to calculate the distance between two points in 3D

def distance(x1, y1, z1, x2, y2, z2):
	dx = x1 - x2
	dy = y1 - y2
	dz = z1 - z2
	dist = math.sqrt(dx*dx + dy*dy + dz*dz)
	return dist 


# Funtion which checks if two atom are close than certain distance and if they are they are added to a neighbouring list

def neigbour(a1, a2, limit, lst):
	x1 = x_co[a1]
        y1 = y_co[a1]
        z1 = z_co[a1]
        x2 = x_co[a2]
        y2 = y_co[a2]
        z2 = z_co[a2]

	distij = distance(x1, y1, z1, x2, y2, z2)
	if distij < limit:
		lst[a1].append(a2)
	return


# Function that finds cycles in the graph
# https://stackoverflow.com/questions/49902879/how-can-i-print-out-the-generator-object-without-using-the-shell

def dfs(graph, start, end):
    fringe = [(start, [])]
    while fringe:
        state, path = fringe.pop()
        if path and state == end:
            yield path
            continue
        for next_state in graph[state]:
            if next_state in path:
                continue
            fringe.append((next_state, path+[next_state]))


# Function that compares lists if they are the same or not

def compare(list1, list2):
    if len(list1) != len(list2):
        return False
    for element in list1:
        if element not in list2:
            return False
    return True

#
# Main program
#

# Recieve the XYZ coordinates

filename = sys.argv[1]
projname = filename.split(".")[0]


#
# Setting the elements and coordinates for later
#

elem = []
x_co = []
y_co = []
z_co = []

with open(filename, "r") as f:
    for line in f:
        data = line.split()
        elem.append(data[0])
        x_co.append(float(data[1]))
        y_co.append(float(data[2]))
        z_co.append(float(data[3]))


#
# Figuring out which atom is which atoms neighbour
#

connect = [[] for _ in range(len(elem))]
lim = 1.65

for i in range(len(elem)):
	for j in range(len(elem)):
        	if elem[i] == "P" or elem[j] == "P" or elem[i] == "S" or elem[j] == "S":
            		lim = 1.82
        	else:
            		lim = 1.65
		if i != j:
			neigbour(i, j, lim, connect)

# Test
# print(connect)

#
# Figuring out the cycles that exist in the molecule
#

graph = {}

for i, neighbors in enumerate(connect):
    graph[i] = set(neighbors)

cycles = [[node]+path  for node in graph for path in dfs(connect, node, node)]

cycles = [l[1:] for l in cycles]

cycles = [l for l in cycles if len(l) != 2]

for i in range(len(cycles)):
	j = 0
	while j < i:
		if compare(cycles[i], cycles[j]):
			cycles.pop(i)
			cycles.insert(i, [])
		j = j + 1

# Test		
# print(cycles)

cycles = [l for l in cycles if len(l) != 0]

#print(cycles)

#
# Start of HOMA/HOMER calculation.
#

with open('HOMA_HOMER_for_excel', 'a') as f: # Making an easy to use excel file
  f.write("HOMA and HOMER values for {}\n".format(projname))
  f.write("Cycle\tHOMA\tHOMER\n")

# Starting the HOMER calculation for each cycle in the molecule.

bond_error_homer = [False]*len(cycles)
atom_error_homer = [False]*len(cycles)
homer = [1]*len(cycles)

for a in range(len(cycles)): 
  cycle = cycles[a] # Definition of  the cycle that we are looking at.

  r_cc_homer = 0 # Setting up the difference sums.
  r_cn_homer = 0
  r_nn_homer = 0
  r_co_homer = 0

  for b in range(len(cycle)):
    if b == len(cycle) - 1:
      c = 0
    else:
      c = b + 1
        
    atom1 = cycle[b] # Figuring out which 2 atom have to be analized.
    atom2 = cycle[c]
		
    x1 = x_co[atom1] # Reading in the coordinates for them.
    x2 = x_co[atom2]
    y1 = y_co[atom1]
    y2 = y_co[atom2]
    z1 = z_co[atom1]
    z2 = z_co[atom2]
		
    dist12 = distance(x1, y1, z1, x2, y2, z2) 

    if elem[atom1] == "C": # Checking what kinf of bond we have.
      if elem[atom2] == "C":
        distdif_homer = dist12 - opt_cc_homer 
        r_cc_homer = r_cc_homer + distdif_homer*distdif_homer
      elif elem[atom2] == "N":
        distdif_homer = dist12 - opt_cn_homer
        r_cn_homer = r_cn_homer + distdif_homer*distdif_homer
      elif elem[atom2] == "O":
        distdif_homer = dist12 - opt_co_homer
        r_co_homer = r_co_homer + distdif_homer*distdif_homer
      else:
        bond_error_homer[a] = True
    elif elem[atom1] == "N":
      if elem[atom2] == "C":
        distdif_homer = dist12 - opt_cn_homer
        r_cn_homer = r_cn_homer + distdif_homer*distdif_homer
      elif elem[atom2] == "N":
        distdif_homer = dist12 - opt_nn_homer
        r_nn_homer = r_nn_homer + distdif_homer*distdif_homer
      else:
        bond_error_homer[a] = True
    elif elem[atom1] == "O":
      if elem[atom2] == "C":
        distdif_homer = dist12 - opt_co_homer
        r_co_homer = r_co_homer + distdif_homer*distdif_homer
      else:
        bond_error_homer[a] = True
    else:
      atom_error_homer[a] = True

  r_tot_homer = a_cc_homer*r_cc_homer + a_cn_homer*r_cn_homer + a_nn_homer*r_nn_homer + a_co_homer*r_co_homer # Calculating HOMER value.
  dif_tot_homer = r_tot_homer / len(cycle)
  homer[a] = 1 - dif_tot_homer
    
# Starting the HOMA calculation for each cycle in the molecule.

bond_error_homa = [False]*len(cycles)
atom_error_homa = [False]*len(cycles)
homa = [1]*len(cycles)
    
for a in range(len(cycles)): 
  cycle = cycles[a] # Definition of  the cycle that we are looking at.

  r_cc_homa = 0 # Setting up the difference sums.
  r_cn_homa = 0
  r_co_homa = 0
  r_cb_homa = 0
  r_cs_homa = 0
  r_cp_homa = 0
  r_nn_homa = 0
  r_no_homa = 0
  r_nb_homa = 0

  for b in range(len(cycle)):
    if b == len(cycle) - 1:
      c = 0
    else:
      c = b + 1

    atom1 = cycle[b] # Figuring out which 2 atom have to be analized.
    atom2 = cycle[c]
		
    x1 = x_co[atom1] # Reading in the coordinates for them.
    x2 = x_co[atom2]
    y1 = y_co[atom1]
    y2 = y_co[atom2]
    z1 = z_co[atom1]
    z2 = z_co[atom2]
		
    dist12 = distance(x1, y1, z1, x2, y2, z2) 

    if elem[atom1] == "C": # Checking what kinf of bond we have.
      if elem[atom2] == "C":
        distdif_homa = dist12 - opt_cc_homa 
        r_cc_homa = r_cc_homa + distdif_homa*distdif_homa
      elif elem[atom2] == "N":
        distdif_homa = dist12 - opt_cn_homer
        r_cn_homa = r_cn_homa + distdif_homa*distdif_homa
      elif elem[atom2] == "O":
        distdif_homa = dist12 - opt_co_homa
        r_co_homa = r_co_homa + distdif_homa*distdif_homa
      elif elem[atom2] == "B":
        distdif_homa = dist12 - opt_cb_homa
        r_cb_homa = r_cb_homa + distdif_homa*distdif_homa
      elif elem[atom2] == "S":
        distdif_homa = dist12 - opt_cs_homa
        r_cs_homa = r_cs_homa + distdif_homa*distdif_homa
      elif elem[atom2] == "P":
        distdif_homa = dist12 - opt_cp_homa
        r_cp_homa = r_cp_homa + distdif_homa*distdif_homa
      else:
        bond_error_homa[a] = True
    elif elem[atom1] == "N":
      if elem[atom2] == "C":
        distdif_homa = dist12 - opt_cn_homa
        r_cn_homa = r_cn_homa + distdif_homa*distdif_homa
      elif elem[atom2] == "N":
        distdif_homa = dist12 - opt_nn_homa
        r_nn_homa = r_nn_homa + distdif_homa*distdif_homa
      elif elem[atom2] == "O":
        distdif_homa = dist12 - opt_no_homa
        r_no_homa = r_no_homa + distdif_homa*distdif_homa
      elif elem[atom2] == "B":
        distdif_homa = dist12 - opt_nb_homa
        r_nb_homa = r_nb_homa + distdif_homa*distdif_homa
      else:
        bond_error_homa[a] = True
    elif elem[atom1] == "O":
      if elem[atom2] == "C":
        distdif_homa = dist12 - opt_co_homa
        r_co_homa = r_co_homa + distdif_homa*distdif_homa
      if elem[atom2] == "N":
        distdif_homa = dist12 - opt_no_homa
        r_no_homa = r_no_homa + distdif_homa*distdif_homa
      else:
        bond_error_homa[a] = True
    elif elem[atom1] == "B":
      if elem[atom2] == "C":
        distdif_homa = dist12 - opt_cb_homa
        r_cb_homa = r_cb_homa + distdif_homa*distdif_homa
      if elem[atom2] == "N":
        distdif_homa = dist12 - opt_nb_homa
        r_nb_homa = r_nb_homa + distdif_homa*distdif_homa
      else:
        bond_error_homa[a] = True
    elif elem[atom1] == "S":
      if elem[atom2] == "C":
        distdif_homa = dist12 - opt_cs_homa
        r_cs_homa = r_cs_homa + distdif_homa*distdif_homa
      else:
        bond_error_homa[a] = True
    elif elem[atom1] == "P":
      if elem[atom2] == "C":
        distdif_homa = dist12 - opt_cp_homa
        r_cp_homa = r_cp_homa + distdif_homa*distdif_homa
      else:
        bond_error_homa[a] = True
    else:
      atom_error_homa[a] = True
            
  r_tot_homa = a_cc_homa*r_cc_homa + a_cn_homa*r_cn_homa + a_co_homa*r_co_homa + a_cp_homa*r_cp_homa + a_cs_homa*r_cs_homa + a_cb_homa*r_cb_homa + a_nn_homa*r_nn_homa + a_no_homa*r_no_homa + a_nb_homa*r_nb_homa  # Calculating HOMA value.
  dif_tot_homa = r_tot_homa / len(cycle)
  homa[a] = 1 - dif_tot_homa

# Reporting the calculated HOMER/HOMA value for the cycle.
	
for a in range(len(cycles)):
  cycle = cycles[a] # Definition of  the cycle that we are looking at.

  cycle_str = " in the following cycle: "
  for d in cycle:
    cycle_str += str(d+1) + " "
  cycle_str_2 = ""
  for d in cycle:
    cycle_str_2 += str(d+1) + " "	

  if bond_error_homer[a] or atom_error_homer[a]:
    if bond_error_homa[a] or atom_error_homa[a]:
      with open('HOMA_HOMER_for_excel', 'a') as f:
        f.write(cycle_str_2 + "\n")
        f.write("ERROR\t")
        f.write("ERROR\n")
    else:
      with open('HOMA_HOMER_for_excel', 'a') as f:
        f.write(cycle_str_2 + "\t")
        f.write("{:.3f}\t".format(homa[a]))
        f.write("ERROR\n")
  else:
    with open('HOMA_HOMER_for_excel', 'a') as f:
        f.write(cycle_str_2 + "\t")
        f.write("{:.3f}\t".format(homa[a]))
        f.write("{:.3f}\n".format(homer[a]))


with open('HOMA_HOMER_for_excel', 'a') as f:
	f.write("\n")
