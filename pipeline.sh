#!/bin/bash

#pipeline

#Standard error messages and input instructions
if [[ $# -ne 3 ]]
then
    echo "usage: bash pipeline.sh -i <protein.pdb> -l <ligand.pdb> -o <output_direcory>"
    echo "FLAGS:"
    echo "  -i: input entire protein structure."
    echo "  -l: input ligand structure."    
    echo "  -o: choose an output directory."
    exit
fi 

#Assigns necessary variable to inputs
PDB_file=$1
LIG=$2
output_directory=$3

#Fixes residues, adds Hydrogens and formats the input protein
bash Pocketanneal-master/Format.sh -i $PDB_file -o FormattedProtein -d Pocketanneal-master/database
cp FormattedProtein/output.pdb ./protein.pdb

#Uses fpocket to determine protein pockets
fpocket -f protein.pdb

#Prepares files for AutoDock
python scripts/prepare_ligand4.py -l $2 -o lig.pdbqt
python scripts/prepare_receptor4.py -r $1 -o protein.pdbqt

#Docks the ligand on each pocket and gets binding affinities
count=0
for poc in `ls ./protein_out/pockets/*.pdb`
	do
	echo "Analyzing pocket $poc..."
	cen=`python scripts/centre.py $poc`
	python scripts/prepare_gpf4.py -r protein.pdbqt -l lig.pdbqt -o $count.gpf -p gridcenter="$cen"
	python scripts/prepare_dpf4.py -r protein.pdbqt -l lig.pdbqt -p ga_num_evals=100000 -p ga_pop_size=100 -p ga_run=15 -o $count.dpf
	autogrid4 -p $count.gpf -l $count.glg
	autodock4 -p $count.dpf -l $count.dlg
	count=$((count +1))
	done

#Selects the best pocket for the particular ligand
mkdir PocketScores
cp *.dlg PocketScores
python scripts/summarize_results4.py -d PocketScores
BestPocketNumber=`cat summary_of_results_1.0|head -n 2|tail -n 1|tr '/' ','|tr ',' '\t'|awk '{print $2}'`

#Removes docking parameter files to save space, comment out if you wish to keep them
rm *.map *.fld *.xyz
rm *.glg *.dlg *.gpf *.dpf

#Prepares files for PocketAnneal
cat PocketScores/$BestPocketNumber.dlg|sed -n '/Cluster Rank = 1/,/ENDMDL/p'|grep "^ATOM" > IP_LIG.pdb
python scripts/make_pocket.py protein_out/pockets/pocket$BestPocketNumber\_atm.pdb protein.pdb IP_LIG.pdb > pocket.pdb

#Runs PocketAnneal and 
python Pocketanneal-master/Pocketanneal.py -P protein.pdb -p pocket.pdb -n 10 -a 2 -o Results -d Pocketanneal-master/database
python scripts/change_residues.py $output_directory/Results/BEST_protein.pdb $output_directory/protein.pdb
python scripts/chain.py $output_directory/Results/BEST_protein.pdb A > Results/BEST_protein_A 
i3Drefine/bin/i3Drefine.sh BEST_protein_A.pdb 1
python scripts/rmsd.py 3DRefiner_Output/RESULT/REFINED_1.pdb Results/BEST_protein_A > RMSD.txt

#Deletes unnecessary folders and files you wish to save space, comment out if you wish to keep them
rm -r FormattedProtein
rm -r protein_out
rm protein.pdbqt lig.pdbqt

#Moves everything to output file
mkdir $output_directory
mv Results $output_directory
mv PocketScores $output_directory
mv 3DRefiner_Output $output_directory
mv IP_LIG.pdb pocket.pdb protein.pdb RMSD.txt summary_of_results_1.0 $output_directory