
#include "../headers/fpocket.h" 

/**

## ----- GENERAL INFORMATION
##
## FILE 					fpocket.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			01-04-08
##
## ----- SPECIFICATIONS
##
##	Top function(s) to use for looking for pockets in a given protein.
##	This function will call successively all function necessary to
##	perform pocket detection using voronoi vertices.
##
##	No output is writen, just the list of pockets are returned.
##
## ----- MODIFICATIONS HISTORY
##
##	09-02-09	(v)  Drop tiny pocket step added
##	28-11-08	(v)  Comments UTD
##	01-04-08	(v)  Added comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## ----- TODO or SUGGESTIONS
##

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
	pockets search_pocket
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	This function will call all functions needed for the pocket finding algorith
	and will return the list of pockets found on the protein.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pdb *pdb : The pdb data of the protein to handle.
	@ s_fparams  : Parameters of the algorithm
   -----------------------------------------------------------------------------
   ## RETURN:
	A chained list of pockets found, sorted according to the current critera
	(the default is a scoring function)
   -----------------------------------------------------------------------------
*/
c_lst_pockets* search_pocket(s_pdb *pdb, s_fparams *params)
{
/*
	clock_t b, e ;
	time_t bt, et ;
*/
	c_lst_pockets *pockets = NULL ;

	/* Calculate and read voronoi vertices comming from qhull */
/*
 	fprintf(stdout,"========= fpocket algorithm begins =========\n") ;

 	fprintf(stdout, "> Calculating vertices ...\n");

	bt = time(NULL) ;
*/
	s_lst_vvertice *lvert = load_vvertices(pdb, params->min_apol_neigh, 
												params->asph_min_size, 
												params->asph_max_size) ;
/*
	et = time(NULL) ;
 	fprintf(stdout, "> Vertices successfully calculated in apox. %f sec.\n",
					(float) (et-bt)) ;
*/
	
	if(lvert == NULL) {
		fprintf(stderr, "! Vertice calculation failed!\n");
		return NULL ;
	}
	/* First clustering */
/* 		fprintf(stdout,"> Basic clustering ...\n");

		b = clock() ;
*/
	pockets = clusterPockets(lvert, params);

	if(pockets) {
		pockets->vertices = lvert ;
/*
		e = clock() ;
		fprintf(stdout, "> Clustering OK in %f sec.\n",
						((double)e - b) / CLOCKS_PER_SEC) ;
*/

	/* Clustering refinment */

/*
		b = clock() ;
		fprintf(stdout,"> Cluster refinment steps: \n");
*/
		reIndexPockets(pockets) ;/* Create index and calculate statistics */
		drop_tiny(pockets) ;	 /* Create index and calculate statistics */
		reIndexPockets(pockets) ;/* Create index and calculate statistics */

/*
		fprintf(stdout,"\t* 2nd refinment step -> clustering : based on barycenters...\n");
*/
		refinePockets(pockets, params) ;	/* Refine clustering (rapid) */
		reIndexPockets(pockets) ;

/*
		fprintf(stdout,"\t* 3rd refinment step -> single linkage clusturing...\n");
*/
		pck_ml_clust(pockets, params);	/* Single Linkage Clustering */
		reIndexPockets(pockets) ;

	/* Descriptors calculation */
/*
		fprintf(stdout,"> Calculating descriptors and score...\n");
		b = clock() ;
*/
		set_pockets_descriptors(pockets);
/*
		e = clock() ;
		fprintf(stdout, "> Descriptors found in %f sec.\n", ((double)e - b) / CLOCKS_PER_SEC) ;

		fprintf(stdout,"> 4th refinment step -> dropping small and polar pockets...\n");
*/

	/* Drop small and too polar binding pockets */
		dropSmallNpolarPockets(pockets, params);
		reIndexPockets(pockets) ;
/*
		e = clock() ;
		fprintf(stdout, "> Refinment OK in %f sec.\n", ((double)e - b) / CLOCKS_PER_SEC) ;
*/

	/* Sorting pockets */
		sort_pockets(pockets, M_SCORE_SORT_FUNCT) ;
		/*sort_pockets(pockets, M_NASPH_SORT_FUNCT) ;*/

		reIndexPockets(pockets) ;
/*
		fprintf(stdout,"===== fpocket algorithm ends =====\n");
*/
	}
	
	return pockets ;
}

