
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

#ifndef DH_DESCR
#define DH_DESCR

/* ---------------------------------INCLUDES--------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

#include "voronoi.h"
#include "voronoi_lst.h"
#include "atom.h"
#include "aa.h"
#include "utils.h"

/* --------------------------------STRUCTURES---------------------------------*/

typedef struct s_desc 
{
	float hydrophobicity_score, /* Hydropathie score - for each aa */
              volume_score,         /* Volume score - for each aa */
              volume,               /* Pocket volume */
              prop_polar_atm,       /* Proportion of polar atoms */
              mean_asph_ray,        /* Mean alpha sphere radius */
              masph_sacc,           /* Mean alpha sphere solvent accessibility */
              apolar_asphere_prop,  /* Proportion of apolar alpha spheres */
              mean_loc_hyd_dens,    /* Mean local hydrophobic density (from alpha spheres) */
              as_density,           /* Pocket density, defined as mean distance between alpha spheres*/
              as_max_dst,           /* Maximum distance between two alpha spheres */
              /* The following descriptors are all normalized using observed
                 values among all pocket found by the algorithm. These
                 are not set in descriptor.c, but in pocket.c as we have to check
                 all pocket first to store boundaries of the descriptor to
                 normalize. */

              flex,                  /* Normalized flexibility - based on B factors - ABUSIVE */
              nas_norm,              /* Normalized number of alpha sphere */
              polarity_score_norm,
              mean_loc_hyd_dens_norm,/* Normalized mean local hydrophobic density */
              prop_asapol_norm,      /* Normalized proportion of apolar alphasphere */
              as_density_norm,
              as_max_dst_norm
        ;
	
	int aa_compo[20] ;	/* Absolute amino acid composition */
	int nb_asph,		/* Number of alpha spheres */
            polarity_score,	/* Polarity score (based on amino acids properties ; see aa.c & aa.h) */
            charge_score ;	/* Sum of all net charges at pH = 7 (see aa.c & aa.h) */
        float as_max_r ;

} s_desc ;

/* ------------------------------PROTOTYPES---------------------------------- */

s_desc* allocate_s_desc(void) ;
void reset_desc(s_desc *desc) ;

void set_descriptors(s_atm **tatoms, int natoms, s_vvertice **tvert, int nvert, s_desc *desc) ;

int get_vert_apolar_density(s_vvertice **tvert, int nvert, s_vvertice *vert) ;
void set_atom_based_descriptors(s_atm **atoms, int natoms, s_desc *desc) ;
void set_aa_desc(s_desc *desc, const char *aa_name) ;


#endif
