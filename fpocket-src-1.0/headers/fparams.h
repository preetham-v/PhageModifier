

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

#ifndef DH_FPARAMS
#define DH_FPARAMS

/* ------------------------------- INCLUDES ----------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>

#include "utils.h"
#include "memhandler.h"

/* ----------------------------- PUBLIC MACROS ------------------------------ */

/* Options of the pocket finder program */
/* standard parameters */

/* Use min alpha sphere radius of : 3.0 */
#define M_MIN_ASHAPE_SIZE_DEFAULT 3.0

/* Use max alpha sphere radius of : 6.0 */
#define M_MAX_ASHAPE_SIZE_DEFAULT 6.0

/* Use first connection distance (see report) : 2.0 */
#define M_CLUST_MAX_DIST 1.73

/*use second connection distance (see report) : 4.5 */
#define M_REFINE_DIST 4.5

/* At least a proportion of  M_REFINE_MIN_NAPOL_AS apolar alpha spheres in 
 * the pocket 0.0 */
#define M_REFINE_MIN_PROP_APOL_AS 0.0

/* Single linkage clustering connection distance 2.5 */
#define M_SLCLUST_MAX_DIST 2.5

/* Minimum number of common neighbours for single linkage clustering 5 */
#define M_SLCLUST_MIN_NUM_NEIGH 2

/* Number of iterations for the Monte Carlo volume calculation 3000 */
#define M_MC_ITER 3000

/* Precision for "exact" volume integration, set to -1 if not used -1 */
#define M_BASIC_VOL_DIVISION -1

/* Minimum number of alpha spheres for a pocket to be kept */
#define M_MIN_POCK_NB_ASPH 36

/* Minimum number of atoms having a low electronegativity in order to declare 
 * an alpha sphere to be apolar 3 */
#define M__MIN_APOL_NEIGH_DEFAULT 3

/* Parameters flags */
#define M_PAR_PDB_FILE 'f'
#define M_PAR_PDB_LIST 'F'
#define M_PAR_MAX_ASHAPE_SIZE 'M'
#define M_PAR_MIN_ASHAPE_SIZE 'm'
#define M_PAR_MIN_APOL_NEIGH 'A'
#define M_PAR_CLUST_MAX_DIST 'D'
#define M_PAR_SL_MAX_DIST 's'
#define M_PAR_SL_MIN_NUM_NEIGH 'n'
#define M_PAR_MC_ITER 'v'
#define M_PAR_BASIC_VOL_DIVISION 'b'
#define M_PAR_MIN_POCK_NB_ASPH 'i'
#define M_PAR_REFINE_DIST 'r'
#define M_PAR_REFINE_MIN_NAPOL_AS 'p'

#define M_FP_USAGE "\n\
***** USAGE (fpocket) *****\n\
\n\
Pocket finding on a pdb - list of pdb - file(s):             \n\
\t./bin/fpocket -f pdb                                       \n\
\t./bin/fpocket -F pdb_list                                  \n\
\nOPTIONS (find standard parameters in brackets)           \n\n\
\t-m (float)  : Minimum radius of an alpha-sphere.      (3.0)\n\
\t-M (float)  : Maximum radius of an alpha-sphere.      (6.0)\n\
\t-A (int)    : Minimum number of apolar neighbor for        \n\
\t              an a-sphere to be considered as apolar.   (3)\n\
\t-i (int)    : Minimum number of a-sphere per pocket.   (30)\n\
\t-D (float)  : Maximum distance for first clustering        \n\
\t              algorithm.                             (1.73)\n\
\t-s (float)  : Maximum distance for single linkage          \n\
\t              clustering                              (2.5)\n\
\t-n (integer): Minimum number of neighbor close from        \n\
\t              each other (not merged otherwise).        (3)\n\
\t-r (float)  : Maximum distance between two pockets         \n\
\t              barycenter (merged otherwise).          (4.5)\n\
\t-p (float)  : Minimum proportion of apolar sphere in       \n\
\t              a pocket (remove otherwise)             (0.0)\n\
\t-v (integer): Number of Monte-Carlo iteration for the      \n\
\t              calculation of each pocket volume.     (2500)\n\
\t-b (integer): Space approximation for the basic method     \n\
\t              of the volume calculation. Not used by       \n\
\t              default (Monte Carlo approximation is)       \n\
\nSee the manual (man fpocket), or the full documentation for\n\
more information.\n\
***************************\n"

/* --------------------------- PUBLIC STRUCTURES ---------------------------- */
/**
	Structure containing all necessary parameters that can be changed by the user.
	This structure is commun to both programs (validation and pocket finding), 
	even if the pocked finding programm doesn't need some parameters.
*/

typedef struct s_fparams
{
	char pdb_path[M_MAX_PDB_NAME_LEN] ;	/* The pdb file */
	char **pdb_lst ;
	int npdb ;
	
	int min_apol_neigh,		 /* Min number of apolar neighbours for an a-sphere 
								to be an apolar a-sphere */
		sl_clust_min_nneigh, /* Min number of neighbours for single linkage 
								clustering */
		nb_mcv_iter,		 /* Number of iteration for the Monte Carlo volume 
								calculation */
		basic_volume_div,	 /* Box division factor for basic volume calculation */
		min_pock_nb_asph ;	 /* Minimump number of alpha spheres per pocket */

	float clust_max_dist,					/* First clustering distance criteria */
		  refine_min_apolar_asphere_prop,	/* Min proportion of apolar alpha 
											spheres for each pocket */
		  sl_clust_max_dist,	/* Single linkage clusturing distance criteria */
		  refine_clust_dist,	/* Refine clustering distance criteria */
		  asph_min_size,	 	/* Minimum size of alpha spheres to keep */
		  asph_max_size ;		/* Maximum size of alpha spheres to keep */


} s_fparams ;

/* ------------------------------- PROTOTYPES ------------------------------- */

s_fparams* init_def_fparams(void) ;
s_fparams* get_fpocket_args(int nargs, char **args) ;

int parse_clust_max_dist(char *str, s_fparams *p) ;
int parse_sclust_max_dist(char *str, s_fparams *p) ;
int parse_sclust_min_nneigh(char *str, s_fparams *p) ;
int parse_min_apol_neigh(char *str, s_fparams *p) ;
int parse_asph_min_size(char *str, s_fparams *p) ;
int parse_asph_max_size(char *str, s_fparams *p) ;
int parse_mc_niter(char *str, s_fparams *p) ;
int parse_basic_vol_div(char *str, s_fparams *p)  ;
int parse_refine_dist(char *str, s_fparams *p)  ;
int parse_refine_minaap(char *str, s_fparams *p)  ;
int parse_min_pock_nb_asph(char *str, s_fparams *p) ;

int is_fpocket_opt(const char opt) ;

void free_fparams(s_fparams *p) ;
void print_pocket_usage(FILE *f) ;
void print_fparams(s_fparams *p, FILE *f) ;

#endif
