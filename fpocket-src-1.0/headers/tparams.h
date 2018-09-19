
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

#ifndef DH_TPARAMS
#define DH_TPARAMS

/* ------------------------------- INCLUDES ----------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>

#include "fparams.h"
#include "utils.h"
#include "memhandler.h"

// -------------------------------------- PUBLIC MACROS ---------------------------------------------

/* Options of the test program */
#define M_PAR_VALID_INPUT_FILE 'L'
#define M_PAR_LIG_NEIG_DIST 'd'
#define M_PAR_P_STATS_OUT 'o'
#define M_PAR_G_STATS_OUT 'e'
#define M_PAR_KEEP_FP_OUTPUT 'k'


/* Write the statistics output of tpocket to : */
#define M_STATS_OUTP "stats_p.txt"
#define M_STATS_OUTG "stats_g.txt"
#define M_MAX_FILE_NAME_LENGTH 300
/* #define M_PAR_P_STATS_OUT 'o' */
/* in order to get the atom set of the pocket, detect around x A of the ligand*/
#define M_LIG_NEIG_DIST 4.0

#define M_TP_USAGE "\
\n***** USAGE (tpocket) *****\n\
\n\
The program needs as input a file containing at each                 \n\
line a pdb file name (apo + complexe), a ligand code                 \n\
(3 letters), all separeted by a tabulation.                          \n\
The format of each line must therefore be:                           \n\n\
{PATH/}APO.pdb  {PATH/}HOLO.pdb  LIG.                                \n\n\
The ligand code is the resname of the ligand atoms in                \n\
the pdb file of the HOLO form of the protein.                        \n\n\
See the manual for more information.                                 \n\n\
Example of command using default parameters:                         \n\
\t./bin/tpocket -L file_path                                         \n\n\
Options:                                                             \n\
\t-e string  : Write global performance to this file                 \n\
\t             Default name: ./stats_g.txt.           (./stats_g.txt)\n\
\t-o string  : Write pocket detailed statistics to .                 \n\
\t             this file Default name: ./stats_p.txt  (./stats_p.txt)\n\
\t-d float   : Distance criteria for the 2 ways to                   \n\
\t             define the actual pocket               (4.0)          \n\n\
Options specific to fpocket are usable too.\n\
See the manual/documentation for mor information.\n\
***************************\n"

/* --------------------------- PUBLIC STRUCTURES -----------------------------*/


typedef struct s_tparams
{
	char **fapo,
		 **fcomplex,
		 **fligan;
	
	char *p_output;
	char *g_output;
	
	char stats_g[128] ; /* M_STATS_OUTG */
	char stats_p[128] ; /* M_STATS_OUTP */

	float lig_neigh_dist ;
	int nfiles ;
	int keep_fpout ;

/* Parameters for the pocket finder program (also needed for validation program...) */

	s_fparams *fpar ;

} s_tparams ;

/* ------------------------------ PROTOTYPES ---------------------------------*/

s_tparams* init_def_tparams(void) ;
s_tparams* get_tpocket_args(int nargs, char **args) ;
int add_list_data(char *str_list_file, s_tparams *par) ;
int add_prot(char *apo, char *complex, char *ligan, s_tparams *par) ;
int parse_lig_neigh_dist(char *str, s_tparams *p) ;

void free_tparams(s_tparams *p);
void print_test_usage(FILE *f) ;
void print_params(s_tparams *p, FILE *f) ;

#endif
