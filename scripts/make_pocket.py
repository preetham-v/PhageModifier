from __future__ import print_function
import sys
import string


if len(sys.argv)!=4:
	print("Usage: "+sys.argv[0]+"python make_pocket.py <pocket> <pdb> <ligand>")
	sys.exit()

Pocket=sys.argv[1]
PDB_file=sys.argv[2]
LIG_file=sys.argv[3]

Pocket_list=[]
for row in open(sys.argv[1], 'r').readlines():
	Pocket_list.append(row)

PDB_list=[]
for row in open(sys.argv[2], 'r').readlines():
	PDB_list.append(row)

Search_Residues=[]

for row in Pocket_list:
	if row[0:4]=="ATOM":
		if row[17:27] not in Search_Residues:
			Search_Residues.append(row[17:27])

LIG_list=[]
for row in open(sys.argv[3], 'r').readlines():
	list1 = list(row)
	list1[21]="X"
	list1[67:87]="                    "
	list1[77]=list1[13]
	row=''.join(list1)
	LIG_list.append(row)


Final_Pocket=[]
for row in PDB_list:
	if row[17:27] in Search_Residues and row[0:4]=="ATOM":
		Final_Pocket.append(row)

for i in Final_Pocket:
	print(i,end="")

for i in LIG_list:
	list2= list(i)
	list2[0:6]="HETATM"
	i=''.join(list2)
	print(i)
