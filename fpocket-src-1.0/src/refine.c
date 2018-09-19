#include "../headers/refine.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					refine.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			01-04-08
##
## ----- SPECIFICATIONS
##
## This file defins several routines which refine the clustering algorithm
## used by fpocket. In particular, we merge pockets too closed from each other
## we drop small pockets, and we perform reindexation after those operations.
##
## ----- MODIFICATIONS HISTORY
##
##	09-02-09	(v)  Drop tiny pocket routine added
##	28-11-08	(v)  Comments UTD 
##	01-04-08	(v)  Added template for comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## ----- TODO or SUGGESTIONS
##
## (v) Improve and optimize the algorithm.
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
	refinePockets
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Refine algorithm: will merge two pockets whose barycenters are
	close together (distance criteria given in params).
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *pockets : The list of pockets.
	@ s_fparams *params      : Parameters
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void refinePockets(c_lst_pockets *pockets, s_fparams *params)
{
	node_pocket *nextPocket;
	node_pocket *curMobilePocket;

	node_pocket *pcur = NULL ;

	float *pbary = NULL, 
		  *mpbary = NULL, dst ;

	if(pockets) {
		pcur = pockets->first ;
		while(pcur) {
			pbary = pcur->pocket->bary ;
			curMobilePocket = pcur->next ;

			while(curMobilePocket) {
				mpbary = curMobilePocket->pocket->bary ;
				nextPocket = curMobilePocket->next ;
				dst = dist(pbary[0], pbary[1], pbary[2], mpbary[0], mpbary[1], mpbary[2]) ;
				if(dst < params->refine_clust_dist) {
				// Merge pockets if barycentres are close to each other
					mergePockets(pcur, curMobilePocket, pockets);		
				}
				curMobilePocket = nextPocket ;
			}

 			pcur = pcur->next ;
		}
	}
	else {
		fprintf(stderr, "! No pocket to refine! (argument NULL: %p).\n", pockets) ;
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	dropSmallNpolarPockets
	-----------------------------------------------------------------------------
   ## SPECIFICATION:
	Refine algorithm: will remove small pockets (depends on the corresponding
	parameters in params), pockets containing less than NB apolar alpha spheres
	(given in params)..
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *pockets : The list of pockets.
	@ s_fparams *params      : Parameters
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void dropSmallNpolarPockets(c_lst_pockets *pockets, s_fparams *params)
{
	double pasph = 0.0 ;
	node_pocket *npcur = NULL,
				*nextPocket1 = NULL;
	s_pocket *pcur = NULL ;
	
	if(pockets) {
		npcur = pockets->first ;
		while(npcur) {
			pcur = npcur->pocket ;
			nextPocket1 = npcur->next ;
			pasph = (float)((float)pcur->nAlphaApol/(float)pcur->v_lst->n_vertices) ;


			if(pcur->v_lst->n_vertices < (size_t) params->min_pock_nb_asph 
				||  pasph <  (params->refine_min_apolar_asphere_prop)){
			/* If the pocket is too small or has not enough apolar alpha
			 * spheres, drop it */
				dropPocket(pockets, npcur);		
			}

			if(pockets->n_pockets <= 0) fprintf(stderr, "! No Pockets Found while refining\n");
			npcur = nextPocket1 ;
		}
	}
	else {
		fprintf(stderr, "! No pockets to drop from (argument NULL: %p).\n", pockets) ;
	}

}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	drop_tiny
	-----------------------------------------------------------------------------
   ## SPECIFICATION:
	Drop really tiny pockets (< 5 alpha spheres)
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *pockets : The list of pockets.
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void drop_tiny(c_lst_pockets *pockets)
{
	node_pocket *npcur = pockets->first,
				*npnext = NULL ;
	while(npcur) {
		npnext = npcur->next ;

		if(npcur->pocket->v_lst->n_vertices < 2){
		/* If the pocket is really small, drop it */
			dropPocket(pockets, npcur);
		}

		if(pockets->n_pockets <= 0) fprintf(stderr, "! No Pockets Found while refining\n");
		npcur = npnext ;
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	reIndexPockets
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Reindex pockets, after dropping and merging operations on pockets and
	recalculate barycentres in the same loop
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *pockets: The list of pockets.
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void reIndexPockets(c_lst_pockets *pockets)
{
	node_vertice *vcur = NULL ;
	node_pocket *pcur = NULL ;
	s_pocket *pock_cur = NULL ;

	int curPocket = 0,
		n_vert ;

	float posSum[3];

	if(pockets && pockets->n_pockets > 0) {
		pcur = pockets->first ;
		while(pcur) {
			curPocket++ ;						//new index counter
			n_vert = 0 ;

			pock_cur = pcur->pocket ;
			pock_cur->bary[0]=0 ;
			pock_cur->bary[1]=0 ;
			pock_cur->bary[2]=0 ;
				
			posSum[0]=0; posSum[1]=0; posSum[2]=0;
			if(pock_cur->v_lst){
				vcur = pock_cur->v_lst->first;
				
				while(vcur){
					posSum[0] += vcur->vertice->x;
					posSum[1] += vcur->vertice->y;
					posSum[2] += vcur->vertice->z;
					n_vert++;
	
					vcur->vertice->resid = curPocket;	//set new index
					vcur = vcur->next ;
				}
			}
			else {
				fprintf(stderr, "! Empty Pocket...something might be wrong over here ;).\n") ; 
			}
			//set new barycentre
			pock_cur->bary[0] = posSum[0]/(float)n_vert;
			pock_cur->bary[1] = posSum[1]/(float)n_vert;
			pock_cur->bary[2] = posSum[2]/(float)n_vert;
			pcur = pcur->next ;
		}
	}
	else {
		fprintf(stderr, "! No pocket to reindex.\n") ;
	}
}

