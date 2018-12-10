import sys
import string
from math import sqrt

if len(sys.argv)!=3:
	print("Usage: "+sys.argv[0]+"python change_residues.py <protein1 PDB> <protein2 PDB>")
	sys.exit()

protein1=sys.argv[1]
protein2=sys.argv[2]

protein1_list=[]
protein2_list=[]

for row in open(sys.argv[1], 'r').readlines():
	protein1_list.append(row)

for row in open(sys.argv[2], 'r').readlines():
	protein2_list.append(row)	


backbone=["CA"]

all_res_1=[]
all_res_2=[]


for row1 in protein1_list:
	if row1[0:4]=="ATOM" and row1[13:15] in backbone:
		residue_1=row1[17:27]
		all_res_1.append(residue_1)

for row2 in protein2_list:
	if row2[0:4]=="ATOM" and row2[13:15] in backbone:
		residue_2=row2[17:27]
		all_res_2.append(residue_2)

for i in range(0,len(all_res_1)-1):
	if all_res_1[i]!=all_res_2[i]:
		print all_res_1[i], all_res_2[i][0:3]



