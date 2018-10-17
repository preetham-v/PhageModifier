#!/bin/bash

#pipeline

#Standard error messages and input instructions
if [[ $# -ne 6 ]]
then
    echo "usage: bash pipeline.sh -i <protein.pdb> -l <ligand.pdb> -o <output_direcory>"
    echo "FLAGS:"
    echo "  -i: input entire protein structure."
    echo "  -l: input ligand structure."    
    echo "  -o: choose an output directory."
    exit
fi 

#definitions:
while getopts ":i::o::l:" opt; do
    case $opt in
    i)
        PDB_file=$OPTARG
        ;;
    l)
        ligand=$OPTARG
        ;;
    o)
        output_directory=$OPTARG
        ;;
    \?)
        echo "error: flag -$OPTARG does not exist." >&2
        exit
        ;;
    :)
        echo "error: flag -$OPTARG requires an argument." >&2
        exit
        ;;
    esac
done

#Fixes residues, adds Hydrogens and formats the input protein
bash Pocketanneal-master/Format.sh -i $PDB_file -o FormattedProtein -d Pocketanneal-master/database
cp FormattedProtein/output.pdb ./protein.pdb

#Uses fpocket to determine protein pockets
fpocket -f protein.pdb

#Prepares files for AutoDock
python scripts/prepare_ligand4.py -l $ligand -o lig.pdbqt
python scripts/prepare_receptor4.py -r protein.pdb -o protein.pdbqt

#Docks the ligand on each pocket and gets binding affinities
count=0
for poc in `ls -v ./protein_out/pockets/*.pdb`
	do
	echo "Analyzing pocket $poc..."
	cen=`python scripts/centre.py $poc`
	python scripts/prepare_gpf4.py -r protein.pdbqt -l lig.pdbqt -o $count.gpf -p gridcenter="$cen"
	python scripts/prepare_dpf4.py -r protein.pdbqt -l lig.pdbqt -p ga_num_evals=100000 -p ga_pop_size=100 -p ga_run=20 -o $count.dpf
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

#Runs PocketAnneal and i3Drefine
python Pocketanneal-master/Pocketanneal.py -P protein.pdb -p pocket.pdb -n 1000 -a 2 -o Results -d Pocketanneal-master/database
python scripts/change_residues.py protein.pdb Results/BEST_protein.pdb > ListOfMutations.txt
python scripts/chain.py Results/BEST_protein.pdb A > Results/BEST_protein_A.pdb 
i3Drefine/bin/i3Drefine.sh Results/BEST_protein_A.pdb 1
python scripts/rmsd.py 3DRefiner_Output/RESULT/REFINED_1.pdb Results/BEST_protein_A.pdb > RMSD.txt

#Deletes unnecessary folders and files to save space, comment out if you wish to keep them
rm -r FormattedProtein
rm -r protein_out
rm protein.pdbqt lig.pdbqt IP_LIG.pdb

#Moves everything to output directory
mkdir $output_directory
mv Results $output_directory
mv PocketScores $output_directory
mv 3DRefiner_Output $output_directory
mv summary_of_results_1.0 $output_directory/Summary_of_pocket_scores
mv pocket.pdb protein.pdb RMSD.txt ListOfMutations.txt $output_directory

echo "Change is the only constant. Till the next protein modification!"