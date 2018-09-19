
#include "../headers/fpout.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					fpout.h
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			28-11-08
##
## ----- SPECIFICATIONS
##
##	Write output for fpocket.
##
## ----- MODIFICATIONS HISTORY
##
##	12-02-09	(v)  No more pocket.info output (useless...)
##	15-12-08	(v)  Minor bug corrected (output dir in the current dir...)
##	28-11-08	(v)  Last argument of write_out_fpocket changed to char *
##					 Comments UTD
##	01-04-08	(v)  Added comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## ----- TODO or SUGGESTIONS
##
##	(v) Handle system command failure, clean!

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
	write_out_fpocket
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Output routine. See the documentation for more information.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
 *  @ c_lst_pockets *pockets : All pockets found and kept.
 *  @ c_lst_pockets *pockets : The (input) pdb structure
	@ char *pdbname          : Name of the pdb
   -----------------------------------------------------------------------------
   ## RETURN: 
	void
   -----------------------------------------------------------------------------
*/
void write_out_fpocket(c_lst_pockets *pockets, s_pdb *pdb, char *pdbname) 
{
	char pdb_code[350] = "" ;
	char pdb_path[350] = "" ;
	char out_path[350] = "" ;
	char pdb_out_path[350] = "" ;
	char fout[350] = "" ;
	char command[370] = "" ;
	int status ;
	
	if(pockets) {
	/* Extract path, pdb code... */
		strcpy(pdb_code, pdbname) ;
		extract_path(pdbname, pdb_path) ;
		remove_ext(pdb_code) ;
		remove_path(pdb_code) ;

		if(strlen(pdb_path) > 0) sprintf(out_path, "%s/%s_out", pdb_path, pdb_code) ;
		else sprintf(out_path, "%s_out", pdb_code) ;
		
		sprintf(command, "mkdir %s", out_path) ;
		status = system(command) ;
		if(status != 0) {
			return ;
		}
		
		sprintf(out_path, "%s/%s", out_path, pdb_code) ;
		sprintf(pdb_out_path, "%s_out.pdb", out_path) ;
		
	/* Write vmd and pymol scripts */
		sprintf(fout, "%s_out.pdb", pdb_code) ;
		write_visualization(out_path, fout);

	/* Writing full pdb */
		sprintf(pdb_out_path, "%s_out.pdb", out_path) ;

		write_pockets_single_pdb(pdb_out_path, pdb, pockets) ;
	
	/* Writing pockets as a single pqr */
		sprintf(fout, "%s_pockets.pqr", out_path) ;
		write_pockets_single_pqr(fout, pockets) ;

	/* Writing individual pockets pqr */
		if(strlen(pdb_path) > 0) sprintf(out_path, "%s/%s_out", pdb_path, pdb_code) ;
		else sprintf(out_path, "%s_out", pdb_code) ;
		
		sprintf(out_path, "%s/pockets", out_path) ;
		sprintf(command, "mkdir %s", out_path) ;
		status = system(command) ;
		if(status != 0) {
			return ;
		}

		write_each_pocket(out_path, pockets) ;
	}
}
