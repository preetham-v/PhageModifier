
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
 
#ifndef DH_WRITE_VISU
#define DH_WRITE_VISU

/* ------------------------------- INCLUDES ----------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "utils.h"

/* ----------------------------- PUBLIC MACROS ------------------------------ */

/* ---------------------------- PUBLIC STRUCTURES --------------------------- */


/* ------------------------------- PROTOTYPES ------------------------------- */


void write_visualization(char *pdb_name,char *pdb_out_name);
void write_vmd(char *pdb_name,char *pdb_out_name);
void write_pymol(char *pdb_name,char *pdb_out_name);

#endif
 
