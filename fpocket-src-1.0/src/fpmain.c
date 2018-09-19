
#include "../headers/fpmain.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					fpmain.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			01-04-08
##
## ----- SPECIFICATIONS
##
##	Top function to call fpocket routines. Get and check programm parameters,
##	call functions, write output and free memory.
##
## ----- MODIFICATIONS HISTORY
##
##	19-01-09	(v)  Minor modif (print on the same line)
##	28-11-08	(v)  process_pdb added, list of pdb taken into account as input
##					 Comments UTD.
##	27-11-08	(v)  PDB file check moved here instead of fparams
##	01-04-08	(v)  Added comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## ----- TODO or SUGGESTIONS
##

*/

/**
    COPYRIGHT DISCLAIMER

    Vincent Le Guilloux, Peter Schmidtke and Pierre Tuffery, hereby
	disclaim all copyright interest in the program “fpocket” (which
	performs protein cavity detection) written by Vincent Le Guilloux and Peter
	Schmidtke.

    Vincent Le Guilloux  28 November 2008
    Peter Schmidtke      28 November 2008
    Pierre Tuffery       28 November 2008

    GNU GPL

    This file is part of the fpocket package.

    fpocket is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    fpocket is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fpocket.  If not, see <http://www.gnu.org/licenses/>.

**/


/**-----------------------------------------------------------------------------
   ## FUNCTION:
	int main(int argc, char *argv[])
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Main program!
   -----------------------------------------------------------------------------
*/
int main(int argc, char *argv[])
{
	fprintf(stdout, "***** POCKET HUNTING BEGINS ***** \n") ;

	s_fparams *params = get_fpocket_args(argc, argv) ;
	
	/* If parameters parsing is ok */
	if(params) {
		if(params->pdb_lst != NULL) {
		/* Handle a list of pdb */
			int i ;
            for (i = 0 ; i < params->npdb ; i++) {
				
				printf("> Protein %d / %d : %s", i, params->npdb,
												   params->pdb_lst[i]) ;
				if(i == params->npdb - 1) fprintf(stdout, "\n") ;
				else fprintf(stdout, "\r") ;
				fflush(stdout) ;
                process_pdb(params->pdb_lst[i], params) ;
            }
		}
		else {
			if(params->pdb_path == NULL || strlen(params->pdb_path) <= 0) {
				fprintf(stdout, "! Invalid pdb name given.\n");
				print_pocket_usage(stdout) ;
			}
			else {
				process_pdb(params->pdb_path, params) ;
			}
		}
	
		free_fparams(params) ;
	}
	else {
		print_pocket_usage(stdout) ;
	}

	fprintf(stdout, "***** POCKET HUNTING ENDS ***** \n") ;
	free_all() ;

	return 0;
}
/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	process_pdb
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Handle a single pdb: check the pdb name, load data, and launch fpocket if
	the pdb file have been successfully read.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ char *pdbname     : Name of the pdb
	@ s_fparams *params : Parameters of the algorithm. See fparams.c/.h
   -----------------------------------------------------------------------------
   ## RETURN: 
	void
   -----------------------------------------------------------------------------
*/
void process_pdb(char *pdbname, s_fparams *params) 
{
	/* Check the PDB file */
	if(pdbname == NULL) return ;
	
	int len = strlen(pdbname) ;
	if(len >= M_MAX_PDB_NAME_LEN || len <= 0) {
		fprintf(stderr, "! Invalid length for the pdb file name. (Max: %d, Min 1)\n",
				M_MAX_PDB_NAME_LEN) ;
		return ;
	}
	
	/* Try to open it */
	s_pdb *pdb =  rpdb_open(pdbname, NULL, M_DONT_KEEP_LIG) ;
	
	if(pdb) {
		/* Actual reading of pdb data and then calculation */
			rpdb_read(pdb, NULL, M_DONT_KEEP_LIG) ;
			c_lst_pockets *pockets = search_pocket(pdb, params);
			if(pockets) {
				write_out_fpocket(pockets, pdb, pdbname);
				c_lst_pocket_free(pockets) ;
			}
	}
	else fprintf(stderr, "! PDB reading failed!\n");
}
