
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

#ifndef DH_RPDBB
#define DH_RPDBB

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "pertable.h"
#include "utils.h"
#include "memhandler.h"


#define M_PDB_LINE_LEN 80   /* actual record size */
#define M_PDB_BUF_LEN  83    /* size need to buffer + CR, LF, and NUL */

#define M_KEEP_LIG  1
#define M_DONT_KEEP_LIG 0

#define M_PDB_HEADER  1
#define M_PDB_REMARK  2
#define M_PDB_ATOM    3
#define M_PDB_CONECT  4
#define M_PDB_HETATM  5
#define M_PDB_CRYST1  6
#define M_PDB_EOF     7
#define M_PDB_END     8
#define M_PDB_UNKNOWN 9

/*
 * API functions start here
 */

typedef struct s_pdb
{
    FILE *fpdb ;

    s_atm *latoms ;     /* The list of atoms: contains all atoms! */

    s_atm **latoms_p ;  /* List of pointers to latoms elements. */
    s_atm **lhetatm ;	/* List of pointer to heteroatoms in the latoms list. */
    s_atm **latm_lig ;	/* List of pointer to the ligand atoms in the atom list*/

    int natoms,			/* Number of atoms */
            nhetatm,		/* Number of HETATM */
            natm_lig ;		/* Number of ligand atoms */

    float A, B, C, 			/* Side lengths of the unit cell */
          alpha, beta, gamma ;	/* Angle between B and C, A and C, A and C */

    char header[M_PDB_BUF_LEN] ;

} s_pdb ;


/* ------------------------------ PUBLIC FUNCTIONS ---------------------------*/

s_pdb* rpdb_open(char *fpath, const char *ligan, const int keep_lig) ;
void rpdb_read(s_pdb *pdb, const char *ligan, const int keep_lig) ;

void rpdb_extract_atm_resname(char *pdb_line, char *res_name) ;

void guess_element(char *aname, char *element) ;

void rpdb_extract_cryst1(char *rstr, float *alpha, float *beta, float *gamma, 
						 float *a, float *b, float *c) ;
void rpdb_extract_atom_values(char *pdb_line, float *x, float *y, float *z,
							  float *occ, float *beta) ;

void rpdb_extract_pdb_atom( char *pdb_line, char *type, int *atm_id, char *name, 
							char *alt_loc, char *res_name, char *chain, 
							int *res_id, char *insert, 
							float *x, float *y, float *z, float *occ, 
							float *bfactor, char *symbol, int *charge, int *guess_flag) ;

void free_pdb_atoms(s_pdb *pdb) ;

#endif
