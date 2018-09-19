
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
 
#ifndef DH_WRITEPOCKET
#define DH_WRITEPOCKET

/* ------------------------------INCLUDES------------------------------------ */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "voronoi.h"
#include "pocket.h"
#include "writepdb.h"
#include "utils.h"

/* ------------------------- PUBLIC STRUCTURES ------------------------------ */

/* -----------------------------PROTOTYPES----------------------------------- */

void write_pockets_single_pdb(const char out[], s_pdb *pdb, c_lst_pockets *pockets)  ;
void write_pockets_single_pqr(const char out[], c_lst_pockets *pockets);

void write_each_pocket(const char out_path[], c_lst_pockets *pockets) ;
void write_pocket_pdb(const char out[], s_pocket *pocket) ;
void write_pocket_pqr(const char out[], s_pocket *pocket) ;

void write_pdb_atoms(FILE *f, s_atm *atoms, int natoms) ;

#endif
