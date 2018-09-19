
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

#ifndef DH_ATOM
#define DH_ATOM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "utils.h"

// --------------------------------MACROS------------------------------------ */



// -------------------------- PUBLIC STRUCTURES ----------------------------- */

/***
	A structure for the modelisation of an atom
**/
typedef struct s_atm
{
    int sort_x;                 /* Index in the sorted tab by X coord */

    float x, y, z ;		/* Coords */
    char name[5],		/* Atom name */
         type[7],		/* Atom type */
         chain[2],		/* Chain name */
         symbol[3],		/* Chemical symbol of the atom */
         res_name[8];		/* Atom residue name */

    int id,			/* Atom id */
        seen,                   /* Say if we have seen this atom during a neighbor search */
        res_id,			/* Atom residue ID */
        atype,
        charge ;		/* Atom charge */

    /* Optional fields */
    float mass,			// Mass */
          radius,		// Vdw radius */
          electroneg,		/* Electronegativity */
          occupancy,		// Occupancy */
          bfactor ;		// B-factor for christal structures */

    char pdb_insert, 		/* PDB insertion code */
         pdb_aloc;		/* PDB alternate location code */

    int atomic_num ;   		/* Atomic number */

} s_atm ;

/* --------------------------------PROTOTYPES--------------------------------- */

float get_mol_mass(s_atm *latoms, int natoms) ;
float get_mol_mass_ptr(s_atm **latoms, int natoms);
void set_mol_barycenter_ptr(s_atm **latoms, int natoms, float bary[3]) ;
float get_mol_volume_ptr(s_atm **atoms, int natoms, int niter) ;

int is_in_lst_atm(s_atm **lst_atm, int nb_atm, int atm_id) ;	
float atm_corsp(s_atm **al1, int nl1, s_atm **pocket_neigh, int nal2) ; 

void print_atoms(FILE *f, s_atm *atoms, int natoms) ;

#endif
