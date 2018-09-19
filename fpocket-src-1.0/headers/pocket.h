
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

#ifndef DH_POCKET
#define DH_POCKET


/* Maximal distance between two vertices in order to belong to the same 
 * pocket (single linkage clustering) */
#define MAX_CON_DIST 2.5

/* -----------------------------INCLUDES--------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "voronoi_lst.h"
#include "pscoring.h"
#include "fparams.h"
#include "aa.h"

#include "utils.h"
#include "memhandler.h"

/* ---------------------- PUBLIC STRUCTURES ----------------------------------*/

typedef struct s_pocket
{
    s_desc *pdesc ;

    c_lst_vertices *v_lst ;
    float score,	/* Discretize different parts of the score */
	  ovlp,         /* Atomic overlap	*/
	  ovlp2,
	  vol_corresp,	/* Volume overlap */
	  bary[3] ;			/* Barycenter of the pocket */
    
    int rank,				/* Rank of the pocket */
	size,
	nAlphaApol,			/* Number of apolar alpha spheres*/
	nAlphaPol ;			/* Number of polar alpha spheres */

} s_pocket ; 

/* Chained list stuff for pockets */
typedef struct node_pocket 
{
	struct node_pocket *next ;
	struct node_pocket *prev ;
	s_pocket *pocket ;

} node_pocket ;

typedef struct c_lst_pockets 
{
	struct node_pocket *first ;
	struct node_pocket *last ;
	struct node_pocket *current ;
	size_t n_pockets ;

	s_lst_vvertice *vertices ;

} c_lst_pockets ;


/* ------------------------------PROTOTYPES----------------------------------*/


/* CLUSTERING FUNCTIONS */
c_lst_pockets *clusterPockets(s_lst_vvertice *lvvert, s_fparams *params);
int updateIds(s_lst_vvertice *lvvert, int i, int *vNb, int resid,int curPocket,c_lst_pockets *pockets, s_fparams *params);
void addStats(int resid, int size, int **stats,int *lenStats);

/* DESCRIPTOR FUNCTIONS */
void set_pockets_descriptors(c_lst_pockets *pockets) ;
void set_normalized_descriptors(c_lst_pockets *pockets) ;
void set_pockets_bary(c_lst_pockets *pockets) ;
s_atm** get_pocket_contacted_atms(s_pocket *pocket, int *natoms) ;
int count_pocket_contacted_atms(s_pocket *pocket) ;
s_vvertice** get_pocket_pvertices(s_pocket *pocket) ;

float set_pocket_mtvolume(s_pocket *pocket, int niter) ;
float set_pocket_volume(s_pocket *pocket, int discret) ;

/* ALLOCATION AND CHAINED LIST OPERATIONS FUNCTIONS*/

s_pocket* alloc_pocket(void) ;
c_lst_pockets *c_lst_pockets_alloc(void);
node_pocket *node_pocket_alloc(s_pocket *pocket);
void c_lst_pocket_free(c_lst_pockets *lst);

node_pocket *c_lst_pockets_add_first(c_lst_pockets *lst, s_pocket *pocket);
node_pocket *c_lst_pockets_add_last(c_lst_pockets *lst,s_pocket *pocket,int cur_n_apol, int cur_n_pol);
void swap_pockets(c_lst_pockets *pockets, node_pocket *p1, node_pocket *p2) ;
void dropPocket(c_lst_pockets *pockets,node_pocket *pocket);
void mergePockets(node_pocket *pocket,node_pocket *pocket2,c_lst_pockets *pockets);
node_pocket *searchPocket(int resid,c_lst_pockets *pockets);

/* OTHER FUNCTON*/
void reset_pocket(s_pocket *pocket) ;
void print_pocket(FILE *f, s_pocket *pocket);
void print_pockets(FILE *f, c_lst_pockets *pockets) ;
void print_pockets_inv(FILE *f, c_lst_pockets *pockets) ;


#endif
