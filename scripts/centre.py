import sys
import string

if len(sys.argv)!=2:
	print("Usage: "+sys.argv[0]+"python centre.py <pdb>")
	sys.exit()

PDB_file=sys.argv[1]

PDB_list=[]
for row in open(sys.argv[1], 'r').readlines():
	PDB_list.append(row)
all_x=[]
all_y=[]
all_z=[]
count=0
sum_x=0
sum_y=0
sum_z=0
for row in PDB_list:
	if row[0:4]=="ATOM":
		x_coord = row[32:38]
		x_coord = x_coord.lstrip()
		y_coord = row[40:47]
		y_coord = y_coord.lstrip()
		z_coord = row[48:55]
		z_coord = z_coord.lstrip()
		sum_x=sum_x+float(x_coord)
		sum_y=sum_y+float(y_coord)
		sum_z=sum_z+float(z_coord)
		count=count+1
mean_x=sum_x/count
mean_y=sum_y/count
mean_z=sum_z/count
print(str(mean_x)+","+str(mean_y)+","+str(mean_z))
#print numpy.mean(all_x)
