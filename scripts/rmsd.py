import sys
import string
from math import sqrt

if len(sys.argv)!=3:
	print("Usage: "+sys.argv[0]+"python rmsd.py <protein1 PDB> <protein2 PDB>")
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

all_x_1=[]
all_y_1=[]
all_z_1=[]

all_x_2=[]
all_y_2=[]
all_z_2=[]

for row1 in protein1_list:
	if row1[0:4]=="ATOM" and row1[13:15] in backbone:
		x_coord = row1[32:38]
		x_coord = float(x_coord.lstrip())
		all_x_1.append(x_coord)
		y_coord = row1[40:47]
		y_coord = float(y_coord.lstrip())
		all_y_1.append(y_coord)
		z_coord = row1[48:55]
		z_coord = float(z_coord.lstrip())
		all_z_1.append(z_coord)


for row2 in protein2_list:
	if row2[0:4]=="ATOM" and row2[13:15] in backbone:
		x_coord = row2[32:38]
		x_coord = float(x_coord.lstrip())
		all_x_2.append(x_coord)
		y_coord = row2[40:47]
		y_coord = float(y_coord.lstrip())
		all_y_2.append(y_coord)
		z_coord = row2[48:55]
		z_coord = float(z_coord.lstrip())
		all_z_2.append(z_coord)

x_sum=0
y_sum=0
z_sum=0

for i in range(0,len(all_x_1)-1):
	x_sum=x_sum+((all_x_1[i]-all_x_2[i])*(all_x_1[i]-all_x_2[i]))
	y_sum=y_sum+((all_y_1[i]-all_y_2[i])*(all_y_1[i]-all_y_2[i]))
	z_sum=z_sum+((all_z_1[i]-all_z_2[i])*(all_z_1[i]-all_z_2[i]))

rmsd=sqrt((x_sum+y_sum+z_sum)/len(all_x_1))

print rmsd
