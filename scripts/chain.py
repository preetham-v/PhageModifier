import sys
import string

if len(sys.argv)!=3:
	print("Usage: "+sys.argv[0]+"python chain.py <protein PDB> <chain ID>")
	sys.exit()

PDB_file=sys.argv[1]
chain_ID=sys.argv[2]

PDB_list=[]

for row in open(sys.argv[1], 'r').readlines():
	PDB_list.append(row)

for row in PDB_list:
	if row[0:4]=="ATOM" and row[21]==chain_ID:
		print row
#print numpy.mean(all_x)
