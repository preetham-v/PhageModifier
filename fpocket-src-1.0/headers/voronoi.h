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


#ifndef DH_VORONOI
#define DH_VORONOI

// ---------------------------------INCLUDES----------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "rpdb.h"
#include "writepdb.h"
#include "calc.h"
#include "utils.h"
#include "../src/qhull/qvoronoi.h"

#include "memhandler.h"

/* ----------------------------------MACROS--------------------------------- */

#define M_VORONOI_SUCCESS 0
/*alpha sphere type - hydrophobic alpha sphere */
#define M_APOLAR_AS 0
/*alpha sphere type - hydrophilic alpha sphere */
#define M_POLAR_AS 1
/*tolerance for coordinate imprecion during alpha sphere search	 */
#define M_PREC_TOLERANCE 1e-4

#define M_BUFSIZE 1e7
/* --------------------------------STRUCTURES-------------------------------- */

typedef struct s_vvertice 
{
	int resid ;
	int id,
            seen,       /* Say if we have seen this vertice during a neighbor search */
            qhullId,
            type ;	/* 0 if apolar contacts, 1 if polar */

	float ray ;	/* Ray of voronoi vertice */
	float x,	/* X coord */
              y,	/* Y coord */
	      z ;	/* Z coord */
	
	int sort_x;		/* Index in the sorted tab by X coord */
	int apol_neighbours;	/* number of neighbouring apolar alpha spheres */

	int vneigh[4] ;
 	s_atm *neigh[4] ;	/* The theorical 4 contacted atoms */
	
	float bary[3] ;         /* Barycenter of the pocket */

} s_vvertice ;

typedef struct s_lst_vvertice
{
	s_vvertice *vertices ;			/* List of voronoi vertices */
        s_vvertice **pvertices ;
        
	int *h_tr;
	int n_h_tr;
	int *tr,
		nvert,
		qhullSize ;

} s_lst_vvertice ;

/* -----------------------------PROTOTYPES----------------------------------- */

s_lst_vvertice* load_vvertices(s_pdb *pdb, int min_apol_neigh, 
				float ashape_min_size, float ashape_max_size) ;
float testVvertice(float xyz[3], int curNbIdx[4], s_atm *atoms, 
				   float min_asph_size, float max_asph_size, 
				   s_lst_vvertice *lvvert);

void set_barycenter(s_vvertice *v) ;
int is_in_lst_vert(s_vvertice **lst_vert, int nb_vert, int v_id) ;
int is_in_lst_vert_p(s_vvertice **lst_vert, int nb_vert, s_vvertice *vert);

void write_pqr_vert(FILE *f, s_vvertice *v) ;
void write_pdb_vert(FILE *f, s_vvertice *v) ;
float get_verts_volume_ptr(s_vvertice **verts, int nvert, int niter) ;

void print_vvertices(FILE *f, s_lst_vvertice *lvvert) ;
void free_vert_lst(s_lst_vvertice *lvvert) ;

#endif
