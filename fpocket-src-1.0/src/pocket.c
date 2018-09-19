 
#include "../headers/pocket.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					pocket.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			01-04-08
##
## ----- SPECIFICATIONS
##	
##		All functions dealing with pockets: clustering, descriptors, allocation...
##
## ----- MODIFICATIONS HISTORY
##
##  10-03-09    (v)  Added a function that count the number of atoms in a pocket.
##	09-02-09	(v)  Normalized maximum distance between two alpha sphere added
##	29-01-09	(v)  Normalized density and polarity score added
##	21-01-09	(v)  Normalized descriptors calculation moved in a single fct
##					 Some variable renamed (refractored)
##	14-01-09	(v)  Added some normalized descriptors
##  28-11-08	(v)  Scoring and sorting moved out of this file + minor 
##					 restructuration + comments almost UTD
##	01-04-08	(v)  Added template for comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## ----- TODO or SUGGESTIONS
##
## (v) Peter: comment culstering function ;) (updateIDs namely)
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


/**
 ================================================================================
 ================================================================================

	FUNCTIONS ON POCKET (clusturing, volume calculation...)

 ================================================================================
 ================================================================================
*/

					/* -------------------- */
					/* CLUSTERING FUNCTIONS */
					/* -------------------- */

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	clusterPockets
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	This function takes in argument a list of vertice, and perform a first
	clusturing algorithm to merge vertices close from each others. The distance
	criteria is in the params struct.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_lst_vvertice *lvvert : The list of vertices.
	@ s_fparams *params      : Parameters
   -----------------------------------------------------------------------------
   ## RETURN:
	list of pockets!
   -----------------------------------------------------------------------------
*/
c_lst_pockets *clusterPockets(s_lst_vvertice *lvvert, s_fparams *params)
{
	int i = -1,
		j = -1,
		curPocketId=1,
		vNb[4] = {0,0,0,0};

	
	s_vvertice *vertices = lvvert->vertices,
			   *vcur = NULL ;

	c_lst_pockets *pockets = c_lst_pockets_alloc();		
		
	for(i=0;i<lvvert->nvert;i++) {
		vcur = vertices + i ;
		for(j=0;j<4;j++) vNb[j] = vcur->vneigh[j];

		curPocketId=updateIds(lvvert,i,vNb,vcur->resid, curPocketId,pockets, params);
	}

	node_pocket *p = pockets->first ;
	while(p) {
		p->pocket->size = p->pocket->v_lst->n_vertices ;
		p = p->next ;
	}

	if(pockets->n_pockets > 0) return pockets ;
	else {
		my_free(pockets) ;
		return NULL ;
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	updateIds
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Update ids.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_lst_vvertice *lvvert : The list of vertices.
	@ int i					 : 
	@ int *vNb				 : 
	@ int resid				 :  
	@ int curPocket			 :
	@ c_lst_pockets *pockets : The current list of pockets.
	@ s_fparams *params		 : Parameters
   -----------------------------------------------------------------------------
   int
   -----------------------------------------------------------------------------
*/
int updateIds(s_lst_vvertice *lvvert, int i, int *vNb, int resid, int curPocket,
			  c_lst_pockets *pockets, s_fparams *params)
{
/* 	s_pocket *curPocket; */
	int filteredIdx,j,z;
	int groupCreatedFlag=0;
	int cur_n_apol=0;
	int cur_n_pol=0;

	s_vvertice *vertices = lvvert->vertices ;
	s_vvertice *vert = &(vertices[i]),
			   *fvert = NULL ;

	for(j=0;j<4;j++) {
		if(vNb[j] < lvvert->qhullSize && vNb[j] > 0) {
			filteredIdx = lvvert->tr[vNb[j]];
			if(filteredIdx!=-1 && filteredIdx < lvvert->nvert){
				fvert = &(vertices[filteredIdx]) ;
				
				if(dist(vert->x, vert->y, vert->z, fvert->x, fvert->y, fvert->z) <= params->clust_max_dist){
					groupCreatedFlag=1;
					/* Add a new pocket */
					if(resid == -1 && fvert->resid==-1){						
						resid=++curPocket;
						vert->resid=resid;
						fvert->resid=curPocket;	
						/* Create a new pocket */
						s_pocket *pocket = alloc_pocket();				
						pocket->v_lst=c_lst_vertices_alloc();
						/* Add vertices to the pocket */
						c_lst_vertices_add_last(pocket->v_lst, vert);				
						c_lst_vertices_add_last(pocket->v_lst, fvert);

						if(vert->type==M_APOLAR_AS) cur_n_apol++;
						else cur_n_pol++;

						if(fvert->type==M_APOLAR_AS)cur_n_apol++;
						else cur_n_pol++;

						c_lst_pockets_add_last(pockets, pocket,cur_n_apol,cur_n_pol);
					}
					/* Add new vertice to existing pocket */
					else if(resid!=-1 && fvert->resid==-1) 						
					{
						fvert->resid=resid;
						node_pocket *pocket=searchPocket(resid,pockets);
						if(pocket) c_lst_vertices_add_last(pocket->pocket->v_lst, fvert);

						if(fvert->type==M_APOLAR_AS) pocket->pocket->nAlphaApol+=1;
						else pocket->pocket->nAlphaPol+=1;
					}
					else if(resid==-1 && fvert->resid!=-1) {
						vert->resid = fvert->resid;
						resid=fvert->resid;
						node_pocket *pocket = searchPocket(resid,pockets);
						if(pocket) c_lst_vertices_add_last(pocket->pocket->v_lst, vert);
						if(vert->type==M_APOLAR_AS) pocket->pocket->nAlphaApol+=1;
						else pocket->pocket->nAlphaPol+=1;
					}
					else if((resid!=-1 && fvert->resid!=-1) && (resid!=fvert->resid)) {
						node_pocket *pocket=searchPocket(resid,pockets);
 						node_pocket *pocket2=searchPocket(fvert->resid,pockets);
						/* Write content of pocket1 into pocket2 */
						mergePockets(pocket,pocket2,pockets);					

						for(z=0;z<lvvert->nvert;z++){
							/* Merge two clusters -> to optimize */
							if(vertices[z].resid==resid) vertices[z].resid=fvert->resid;	
						}
						resid=fvert->resid;

					}
				}
			}
			else {
				if(filteredIdx >= lvvert->nvert) {
					fprintf(stdout, "INDEX ERROR\n") ;
				}
			}
		}
	}
	/* if a pocket was created return the pocket */
	if(groupCreatedFlag) return curPocket;
	
	/* if no vertice neighbours were found, still build a new pocket containing 
	 * only one vertice */
	
	cur_n_apol=0;
	cur_n_pol=0;
	resid=++curPocket;
	vert->resid=resid;
	s_pocket *pocket=alloc_pocket();				/* Create a new pocket */
	pocket->v_lst=c_lst_vertices_alloc();
	c_lst_vertices_add_last(pocket->v_lst, vert);	/* Add vertices to the pocket */
	if(vert->type==M_APOLAR_AS) cur_n_apol++;
	else cur_n_pol++;

	c_lst_pockets_add_last(pockets,pocket,cur_n_apol,cur_n_pol);

	return curPocket;
	
}


				/* ------------------------------- */
				/* VOLUME AND DESCRIPTOR FUNCTIONS */
				/* ------------------------------- */


/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	set_pocket_mtvolume
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Get an monte carlo approximation of the volume occupied by the pocket given
	in argument.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pocket *pocket: Pockets to print
	@ int niter       : Number of monte carlo iteration
   -----------------------------------------------------------------------------
   ## RETURN:
	float: volume.
   -----------------------------------------------------------------------------
*/
float set_pocket_mtvolume(s_pocket *pocket, int niter) 
{
	int i = 0,
		nb_in = 0,
		nit = 0 ;

	float xmin = 0.0, xmax = 0.0,
		  ymin = 0.0, ymax = 0.0,
		  zmin = 0.0, zmax = 0.0,
		  xtmp = 0.0, ytmp = 0.0, ztmp = 0.0,
		  xr = 0.0, yr = 0.0, zr = 0.0,
		  vbox = 0.0 ;

	c_lst_vertices *vertices = pocket->v_lst ;
	node_vertice *cur = vertices->first ;
	s_vvertice *vcur = NULL ;

	/* First, search extrems coordinates */
	while(cur) {
		vcur = cur->vertice ;
		
		/* Update min: */
		if(nit == 0) {
			xmin = vcur->x - vcur->ray ; xmax = vcur->x + vcur->ray ;
			ymin = vcur->y - vcur->ray ; ymax = vcur->y + vcur->ray ;
			zmin = vcur->z - vcur->ray ; zmax = vcur->z + vcur->ray ;
		}
		else {
			if(xmin > (xtmp = vcur->x - vcur->ray)) xmin = xtmp ;
			else if(xmax < (xtmp = vcur->x + vcur->ray)) xmax = xtmp ;
	
			if(ymin > (ytmp = vcur->y - vcur->ray)) ymin = ytmp ;
			else if(ymax < (ytmp = vcur->y + vcur->ray)) ymax = ytmp ;
	
			if(zmin > (ztmp = vcur->z - vcur->ray)) zmin = ztmp ;
			else if(zmax < (ztmp = vcur->z + vcur->ray)) zmax = ztmp ;
		}

		cur = cur->next ;
		nit++ ;
	}

	/* Next calculate the box volume */
	vbox = (xmax - xmin)*(ymax - ymin)*(zmax - zmin) ;
	
	/* Then apply monte carlo approximation of the volume.	 */
	for(i = 0 ; i < niter ; i++) {
		xr = rand_uniform(xmin, xmax) ;
		yr = rand_uniform(ymin, ymax) ;
		zr = rand_uniform(zmin, zmax) ;
		cur = vertices->first ;
		
		while(cur) {
			vcur = cur->vertice ;
		/* Distance between the center of curent vertice and the random point */
			xtmp = vcur->x - xr ;
			ytmp = vcur->y - yr ;
			ztmp = vcur->z - zr ;

		/* Compare r^2 and dist(center, random_point)^2  to avoid a call to sqrt() function */
			if((vcur->ray*vcur->ray) > (xtmp*xtmp + ytmp*ytmp + ztmp*ztmp)) {
			/* The point is inside one of the vertice!! */
				nb_in ++ ; break ;
			}
			cur = cur->next ;
		}
	}

	pocket->pdesc->volume = ((float)nb_in)/((float)niter)*vbox ;

	/* Ok lets just return the volume Vpok = Nb_in/Niter*Vbox */
	return pocket->pdesc->volume ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	set_pocket_volume
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Get an approximation of the volume occupied by the pocket given in argument,
	using a discretized space.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pocket *pocket: Pockets to print
	@ int idiscret    : Discretisation: the cube containing all vertices will be 
						divided in idiscret*idiscret*idiscret cubes.
   -----------------------------------------------------------------------------
   ## RETURN:
	float: volume.
   -----------------------------------------------------------------------------
*/
float set_pocket_volume(s_pocket *pocket, int idiscret) 
{
	int niter = 0,
		nb_in = 0,
		nit = 0 ;

	float discret = 1.0/(float)idiscret ;

	float x = 0.0, y = 0.0, z = 0.0,
		  xstep = 0.0, ystep = 0.0, zstep = 0.0 ;

	float xmin = 0.0, xmax = 0.0,
		  ymin = 0.0, ymax = 0.0,
		  zmin = 0.0, zmax = 0.0,
		  xtmp = 0.0, ytmp = 0.0, ztmp = 0.0,
		  vbox = 0.0 ;

	c_lst_vertices *vertices = pocket->v_lst ;
	node_vertice *cur = vertices->first ;
	s_vvertice *vcur = NULL ;

	/* First, search extrems coordinates */
	while(cur) {
		vcur = cur->vertice ;
		
		/* Update min: */
		if(nit == 0) {
			xmin = vcur->x - vcur->ray ; xmax = vcur->x + vcur->ray ;
			ymin = vcur->y - vcur->ray ; ymax = vcur->y + vcur->ray ;
			zmin = vcur->z - vcur->ray ; zmax = vcur->z + vcur->ray ;
		}
		else {
			if(xmin > (xtmp = vcur->x - vcur->ray)) xmin = xtmp ;
			else if(xmax < (xtmp = vcur->x + vcur->ray)) xmax = xtmp ;
	
			if(ymin > (ytmp = vcur->y - vcur->ray)) ymin = ytmp ;
			else if(ymax < (ytmp = vcur->y + vcur->ray)) ymax = ytmp ;
	
			if(zmin > (ztmp = vcur->z - vcur->ray)) zmin = ztmp ;
			else if(zmax < (ztmp = vcur->z + vcur->ray)) zmax = ztmp ;
		}

		cur = cur->next ;
		nit++ ;
	}

	/* Next calculate the box volume */
	vbox = (xmax - xmin)*(ymax - ymin)*(zmax - zmin) ;

	xstep = discret * (xmax - xmin) ;
	ystep = discret * (ymax - ymin) ;
	zstep = discret * (zmax - zmin) ;

	/* Then apply monte carlo approximation of the volume.	 */
	for(x = xmin ; x < xmax ; x += xstep) {
		for(y = ymin ; y < ymax ; y += ystep) {	
			for(z = zmin ; z < zmax ; z += zstep) {
				cur = vertices->first ;
				while(cur) {
					vcur = cur->vertice ;
					xtmp = vcur->x - x ;
					ytmp = vcur->y - y ;
					ztmp = vcur->z - z ;
		
				/* Compare r^2 and dist(center, random_point)^2 */
					if((vcur->ray*vcur->ray) > (xtmp*xtmp + ytmp*ytmp + ztmp*ztmp)) {
					/*the point is inside one of the vertice!! */
						nb_in ++ ; break ;
					}
					cur = cur->next ;
				}
				niter ++ ;
			}
		}
	}

	pocket->pdesc->volume = ((float)nb_in)/((float)niter)*vbox ;

	/* Ok lets just return the volume Vpok = Nb_in/Niter*Vbox */
	return pocket->pdesc->volume ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	set_pockets_bary
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Set barycenter of each pockets. Use vertices to do so.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *pockets: The list of pockets to handle
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void set_pockets_bary(c_lst_pockets *pockets)
{
	node_pocket *cur = NULL ;
	s_pocket *pcur = NULL ;

	float xsum, ysum, zsum ;
	int n ;

	if(pockets && pockets->n_pockets > 0) {
		cur = pockets->first ;
		while(cur) {
			pcur = cur->pocket ;
		/* Reset values and calculate barycenter */
			xsum = 0.0 ; ysum = 0.0 ; zsum = 0.0 ; 
			n = 0 ;

			node_vertice *nvcur = pcur->v_lst->first ;
			while(nvcur) {
				xsum += nvcur->vertice->x ;
				ysum += nvcur->vertice->y ;
				zsum += nvcur->vertice->z ;
				n ++ ;

				nvcur = nvcur->next ;
			}

			pcur->bary[0] = xsum / (float) n ;
			pcur->bary[1] = ysum / (float) n ;
			pcur->bary[2] = zsum / (float) n ;

			cur = cur->next ;
		}
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	set_descriptors
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Set descriptors for a list of pockets.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *pockets: The list of pockets to handle
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void set_pockets_descriptors(c_lst_pockets *pockets)
{
	node_pocket *cur = NULL ;
	s_pocket *pcur = NULL ;
	
	int natms = 0, i ;

	if(pockets && pockets->n_pockets > 0) {
		cur = pockets->first ;
		/* Perform a first loop to calculate atom and vertice based descriptors */
		while(cur) {
			pcur = cur->pocket ;

			/* Get a list of vertices in a tab of pointer */
			s_vvertice **tab_vert = (s_vvertice **)
					my_malloc(pcur->v_lst->n_vertices*sizeof(s_vvertice*)) ;
			i = 0 ;
			node_vertice *nvcur = pcur->v_lst->first ;
			while(nvcur) {
				tab_vert[i] = nvcur->vertice ;
				nvcur = nvcur->next ;
				i++ ;
			}

			/* Get atoms contacted by vertices, and calculate descriptors */
			s_atm **pocket_atoms = get_pocket_contacted_atms(pcur, &natms) ;

			/* Calculate descriptors*/
			set_descriptors(pocket_atoms, natms, tab_vert,
							pcur->v_lst->n_vertices, pcur->pdesc) ;

			my_free(pocket_atoms) ;
			my_free(tab_vert) ;

			cur = cur->next ;
		}

		/* Set normalized descriptors */
		set_normalized_descriptors(pockets) ;

		/* Score all pockets */
		cur = pockets->first ;
		while(cur) {
			cur->pocket->score = score_pocket(cur->pocket->pdesc) ;

			cur = cur->next ;
		}
	}
}


/**-----------------------------------------------------------------------------
   ## FUNCTION:
	set_normalized_descriptors
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	 Perform normalisation for some descriptors, so that the maximum value of a
	 given descriptor become 1 and the minimum value 0. This way, each descriptor
	 is normalized and take into account relative differences between pockets.
	 To do so, we use the simple formula:
		Dnorm = (D - Dmin) / (Dmax - Dmin)
	 with Dmin (resp. dmax) being the minimum (resp. maximum) value of the
 	 descriptor D observed in all detected pockets.

     WARNING: It is assumed that basic descriptors have been normalized before
	 calling this function!
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *pockets : The list of pockets
   -----------------------------------------------------------------------------
   ## RETURN:
    void: s_desc is filled
   -----------------------------------------------------------------------------
*/
void set_normalized_descriptors(c_lst_pockets *pockets)
{
	
	if(!pockets || pockets->n_pockets <= 0) return ;

	node_pocket *cur = NULL ;
	s_desc *dcur = NULL ;

	/* Some boundaries to help normalisation */
	float flex_M = 0.0, flex_m = 1.0,
		  nas_apolp_M = 0.0, nas_apolp_m = 1.0,
		  density_M = 0.0, density_m = 100.0,
		  mlhd_M = 0.0, mlhd_m = 1000.0,
		  as_max_dst_M = 0.0, as_max_dst_m = 1000.0 ;

	int nas_M = 0, nas_m = 1000,
		polarity_M = -1, polarity_m = 10000 ;

	cur = pockets->first ;
	/* Perform a first processing step, to set min and max for example */
	while(cur) {
		dcur = cur->pocket->pdesc ;
		/* Initialize boundaries if it's the first turn */
		if(cur == pockets->first) {
			as_max_dst_M = as_max_dst_m = dcur->as_max_dst ;
			density_M = density_m = dcur->as_density ;
			polarity_M = polarity_m = dcur->polarity_score ;
			flex_M = flex_m = dcur->flex ;
			nas_apolp_M = nas_apolp_m = dcur->apolar_asphere_prop ;
			mlhd_M = mlhd_m = dcur->mean_loc_hyd_dens ;
			nas_M = nas_m = dcur->nb_asph ;
		}
		else {
		/* Update several boundaries */
			if(dcur->as_max_dst > as_max_dst_M)
				as_max_dst_M = dcur->as_max_dst ;
			else if(dcur->as_max_dst < as_max_dst_m)
				as_max_dst_m = dcur->as_max_dst ;

			if(dcur->as_density > density_M)
				density_M = dcur->as_density ;
			else if(dcur->as_density < density_m)
				density_m = dcur->as_density ;

			if(dcur->polarity_score > polarity_M)
				polarity_M = dcur->polarity_score ;
			else if(dcur->polarity_score < polarity_m)
				polarity_m = dcur->polarity_score ;

			if(dcur->mean_loc_hyd_dens > mlhd_M)
				mlhd_M =dcur->mean_loc_hyd_dens ;
			else if(dcur->mean_loc_hyd_dens < mlhd_m)
				mlhd_m =dcur->mean_loc_hyd_dens ;

			if(dcur->flex > flex_M) flex_M = dcur->flex ;
			else if(dcur->flex < flex_m) flex_m = dcur->flex ;

			if(dcur->nb_asph > nas_M) nas_M = dcur->nb_asph ;
			else if(dcur->nb_asph < nas_m) nas_m = dcur->nb_asph ;

			if(dcur->apolar_asphere_prop > nas_apolp_M)
				nas_apolp_M = dcur->apolar_asphere_prop ;
			else if(dcur->apolar_asphere_prop < nas_apolp_m)
				nas_apolp_m = dcur->apolar_asphere_prop;
		}

		cur = cur->next ;
	}

	/* Perform a second loop to do the actual normalisation */
	cur = pockets->first ;
	while(cur) {
		dcur = cur->pocket->pdesc ;
		/* Calculate normalized descriptors */
		dcur->as_max_dst_norm =
			(dcur->as_max_dst - as_max_dst_m) / (as_max_dst_M - as_max_dst_m) ;

		dcur->as_density_norm =
			(dcur->as_density - density_m) / (density_M - density_m) ;

		dcur->polarity_score_norm =
			(float)(dcur->polarity_score - polarity_m) /
			(float)(polarity_M - polarity_m) ;

		dcur->mean_loc_hyd_dens_norm =
			(dcur->mean_loc_hyd_dens - mlhd_m) / (mlhd_M - mlhd_m) ;

		dcur->flex = (dcur->flex - flex_m) / (flex_M - flex_m) ;
		dcur->nas_norm = (float) (dcur->nb_asph - nas_m) /
								(float) (nas_M - nas_m) ;
		dcur->prop_asapol_norm =
			(dcur->apolar_asphere_prop - nas_apolp_m)
			/ (nas_apolp_M - nas_apolp_m);

		cur = cur->next ;
	}
}


/**
 ================================================================================
 ================================================================================

	FUNCTIONS FOR POCKET CHAINED LIST OPERATIONS

 ================================================================================
 ================================================================================
*/

/**-----------------------------------------------------------------------------
   ## FONCTION: 
  	c_lst_pockets_add_first
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Add a pocket on the first position of the list
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pocket *lst : chained list of pockets
	@ s_pocket			: pocket to add
   -----------------------------------------------------------------------------
   ## RETURN:
	node_pocket *: pointer to the new node.
   -----------------------------------------------------------------------------
*/
node_pocket *c_lst_pockets_add_first(c_lst_pockets *lst, s_pocket *pocket)
{
	node_pocket *newn = NULL ;

	if(lst) {
		newn = node_pocket_alloc(pocket) ;
		lst->first->prev = newn ;
		newn->next = lst->first ;

		lst->first = newn ;
		lst->n_pockets += 1 ;
	}
	
	return newn ;
}

/**-----------------------------------------------------------------------------
   ## FONCTION: 
	c_lst_pockets_add_last
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Add a pocket at the end of the chained list
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pocket *lst : chained list of pockets
	@ s_pocket          : pocket to add
   -----------------------------------------------------------------------------
   ## RETURN:
	node_pocket *: Pointer to the new pocket.
   -----------------------------------------------------------------------------
*/
node_pocket *c_lst_pockets_add_last(c_lst_pockets *lst,s_pocket *pocket,int cur_n_apol, int cur_n_pol)
{
	node_pocket *newn = NULL ;

	if(lst) {
		newn = node_pocket_alloc(pocket) ;
		newn->pocket->nAlphaApol = cur_n_apol;
 		newn->pocket->nAlphaPol = cur_n_pol;
		if(lst->last) {
			newn->prev = lst->last ;
			lst->last->next = newn ;
		}
		else {
			lst->first = newn ;

		}
		lst->last = newn ;
		lst->n_pockets += 1 ;
	}

	return newn ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	swap_pockets 
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Swap two pockets in the given list (not that easy !!! ;) )
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *pockets : The list of pockets
	@ const node_pocket *p1  : Pocket 1
	@ const node_pocket *p2  : Pocket 2
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void swap_pockets(c_lst_pockets *pockets, node_pocket *p1, node_pocket *p2) 
{
	node_pocket *p1prev = p1->prev,
				*p1next = p1->next,
				*p2prev = p2->prev,
				*p2next = p2->next ;
	
	if(p1->next == p2) {
	/* If p1 is just before p2 (p1 can't be the last one and p2 can't be
	 * the first one) */
		if(p1->prev) p1->prev->next = p2 ;
		else pockets->first = p2 ; /* P1 is the first of the list */

		if(p2->next) p2->next->prev = p1 ;
		else pockets->last = p2 ; /* P2 is the last one */

		p2->next = p1 ;
		p2->prev = p1prev ;
		p1->prev = p2 ;
		p1->next = p2next ;
	}
	else if(p1->prev == p2) {
	/* If p1 is just after p2 (p2 can't be the last one and p1 can't be
	 * the first one) */
		if(p1->next) p1->next->prev = p2 ;
		else pockets->last = p2 ; /* p1 is the last one */
		
		if(p2->prev) p2->prev->next = p1 ;
		else pockets->first = p1 ;

		p1->next = p2 ;
		p1->prev = p2prev ;
		p2->prev = p1 ;
		p2->next = p1next;
	}
	else {
	/* If p1 and p2 are far away from each others in the list */
	/* Handle p1 */
		if(p1->next) {
		/* If there is a pocket after p1, p1 is not the last one */
			p1->next->prev = p2 ;	/* p2 should be the new prev of p1's next */
			if(p1->prev) p1->prev->next = p2 ;
			else pockets->first = p2 ;
		}
		else {
		/* p1 is the last pocket of the list */
			
			pockets->last = p2 ;
			if(p1->prev) {
				p1->prev->next = p2 ;
			}
		}

	/* Handle P2 */
		if(p2->next) {
		/* If there is a pocket after p1, p1 is not the last one */
			p2->next->prev = p1 ;	/* p2 should be the new prev of p1's next */
			if(p2->prev) p2->prev->next = p1 ;
			else pockets->first = p1 ;
		}
		else {
		/* p1 is the last pocket of the list */
			if(p2->prev) {
				p2->prev->next = p1 ;
				pockets->last = p1 ;
			}
		}

		p1->next = p2next ;
		p1->prev = p2prev ;

		p2->next = p1next ;
		p2->prev = p1prev ;
		
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	searchPocket
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Search a pocket
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ int resid			 : ID of the pocket
	@ c_lst_pockets *lst : The list of pockets
   -----------------------------------------------------------------------------
   ## RETURN:
	node_pocket *: pointer to the pocket found. 
   -----------------------------------------------------------------------------
*/
node_pocket *searchPocket(int resid, c_lst_pockets *lst)
{
	node_pocket *cur = NULL ;
	
	if(lst) {
		cur = lst->first ;
		while(cur && cur->pocket->v_lst->first->vertice->resid != resid) {
			cur = cur->next ;
        }
	}

	return cur;
}

/**-----------------------------------------------------------------------------
   ## FONCTION: 
	void dropPocket(c_lst_pockets *pockets,node_pocket *pocket)
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Drop a pocket from the list.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *lst  : The list.
	@ node_pocket *pocket : The pocket to drop
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void dropPocket(c_lst_pockets *pockets, node_pocket *pocket)
{
	
	if(pocket->prev && pocket->next){
		pocket->prev->next = pocket->next;
		pocket->next->prev = pocket->prev;
	}
	else if(pocket->next){
		pocket->next->prev = NULL;
		pockets->first = pocket->next;
	}
	else if(pocket->prev){
		pocket->prev->next = NULL;
		pockets->last = pocket->prev;
	}

	pocket->next = NULL ;
	pocket->prev = NULL ;
	pockets->n_pockets -= 1;

	my_free(pocket->pocket->pdesc);
	pocket->pocket->pdesc = NULL ;
	my_free(pocket->pocket);
	pocket->pocket= NULL;

	my_free(pocket);

	if(pockets->n_pockets == 0) pockets->first = NULL ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	 mergePockets
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Merge two pockets.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ node_pocket *pocket: The first pocket
	@ node_pocket *pocket2: The second pocket
	@ c_lst_pockets *pockets: The list of pockets
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void mergePockets(node_pocket *pocket,node_pocket *pocket2,c_lst_pockets *pockets)
{
	s_pocket *pock =  pocket->pocket,
			 *pock2 = pocket2->pocket ;

	pock->nAlphaApol += pock2->nAlphaApol;
	pock->nAlphaPol += pock2->nAlphaPol;
	pock->v_lst->n_vertices += pock2->v_lst->n_vertices;
	pock->size = pock->v_lst->n_vertices ;
	
	pock->v_lst->last->next = pock2->v_lst->first;
	pock->v_lst->last = pock2->v_lst->last;
	pocket2->pocket->v_lst = NULL;

	my_free(pocket2->pocket);
	pocket2->pocket=NULL;
	
	if(pocket2->prev && pocket2->next){
		pocket2->prev->next = pocket2->next;
		pocket2->next->prev = pocket2->prev;
	}
	else if(pocket2->next){
		pockets->first=pocket2->next;
		pocket2->next->prev=NULL;
	}
	else if(pocket2->prev){

		pocket2->prev->next=NULL;
		pockets->last=pocket2->prev;
	}
	else if(!pocket2->prev && !pocket2->next) my_free(pocket2);

	pockets->n_pockets-=1;
 	my_free(pocket2);
}

/**
 ================================================================================
 ================================================================================

	ALLOCATIONS AND DESALLOCATIONS FOR POCKET AND CHAINED LIST OF POCKET 

 ================================================================================
 ================================================================================
*/

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	alloc_pocket 
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Alloc memory for a pocket and reste values describing it.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
   -----------------------------------------------------------------------------
   ## RETURN:
	s_pocket* : pointer to the pocket allocated. 
   -----------------------------------------------------------------------------
*/
s_pocket* alloc_pocket(void) 
{
	s_pocket *p = (s_pocket*)my_malloc(sizeof(s_pocket)) ;
	p->pdesc = (s_desc*)my_malloc(sizeof(s_desc)) ;
	p->v_lst = NULL ;

	reset_pocket(p) ;

	return p ;
}

/**-----------------------------------------------------------------------------
   ## FONCTION: 
   c_lst_pockets_alloc
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Allocate a list of pockets
   -----------------------------------------------------------------------------
   ## PARAMETRES:
   -----------------------------------------------------------------------------
   ## RETURN:
	c_lst_pockets*
   -----------------------------------------------------------------------------
*/
c_lst_pockets *c_lst_pockets_alloc(void) 
{
	c_lst_pockets *lst = (c_lst_pockets *)my_malloc(sizeof(c_lst_pockets)) ;

	lst->first = NULL ;
	lst->last = NULL ;
	lst->current = NULL ;
	lst->n_pockets = 0 ;
	lst->vertices = NULL ;

	return lst ;
}

/**-----------------------------------------------------------------------------
   ## FONCTION: 
	node_pocket_alloc
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Allocate memory for one pocket node
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pocket *pocket : pointer to the pocket
   -----------------------------------------------------------------------------
   ## RETURN:
	node_pocket*
   -----------------------------------------------------------------------------
*/
node_pocket *node_pocket_alloc(s_pocket *pocket)
{
	node_pocket *n_pocket = (node_pocket*)my_malloc(sizeof(node_pocket)) ;

	n_pocket->next = NULL ;
	n_pocket->prev = NULL ;
	n_pocket->pocket = pocket ;
	
	return n_pocket ;
}


/**-----------------------------------------------------------------------------
   ## FONCTION: 
	c_lst_pocket_free
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Free a pocket list.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ c_lst_pockets *lst: The list to free.
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void c_lst_pocket_free(c_lst_pockets *lst) 
{
	node_pocket *cur = NULL,
				*next = NULL ;
	
	if(lst) {
		cur = lst->first ;
		while(cur) {
			next = cur->next ;
			dropPocket(lst, cur) ;
			cur = next ;
        }

		free_vert_lst(lst->vertices) ;

		lst->vertices = NULL ;
		lst->first = NULL ;	
		lst->last = NULL ;
		lst->current = NULL ;

		my_free(lst) ;
	}

}

				// ----------------------------------------------- //
				// SORTING FUNCTION (REARANGE CHAINED LIST TO HAVE //
				// POCKETS SORTED ACCORDING TO SEVERAL PROPERTIES) //
				// ----------------------------------------------- //

/**
 ================================================================================
 ================================================================================

	OTHER FUNCTIONS (contacted atoms, print...)

 ================================================================================
 ================================================================================
*/

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	get_pocket_contacted_atms
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Get all atoms contacted by the alpha spheres of a given pocket.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pocket *pocket : The pocket
	@ int *natoms      : OUTPUT Number of atoms found (modified)
   -----------------------------------------------------------------------------
   ## RETURN:
	s_atm** Modify the value of natoms, and return an array of pointer to atoms.
   -----------------------------------------------------------------------------
*/
s_atm** get_pocket_contacted_atms(s_pocket *pocket, int *natoms)
{
	int actual_size = 10,
		nb_atoms = 0,
		i = 0 ;
	
	node_vertice *nvcur = NULL ;
	s_vvertice *vcur = NULL ;

	s_atm **catoms = NULL ;
	
	if(pocket && pocket->v_lst && pocket->v_lst->n_vertices > 0) {
	/* Remember atoms already stored. */
		int atm_ids[pocket->v_lst->n_vertices * 4] ;

	/* Do the search  */
		catoms = (s_atm **)my_malloc(actual_size*sizeof(s_atm*)) ;
		nvcur = pocket->v_lst->first ;
		while(nvcur) {
			vcur = nvcur->vertice ;
			/*printf("ID in the pocket: %d (%.3f %.3f %.3f\n", vcur->id, vcur->x, vcur->y, vcur->z) ;*/
			for(i = 0 ; i < 4 ; i++) {
				if(in_tab(atm_ids,  nb_atoms, vcur->neigh[i]->id) == 0) {
					if(nb_atoms >= actual_size) {
						catoms = (s_atm **)my_realloc(catoms, (actual_size+10)*sizeof(s_atm**)) ;
						actual_size += 10 ;
					}
	
					atm_ids[nb_atoms] = vcur->neigh[i]->id ;
					catoms[nb_atoms] = vcur->neigh[i] ;
					nb_atoms ++ ;
				}
			}
			nvcur = nvcur->next ;
		}
	}

	*natoms = nb_atoms ;
	
	return catoms ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	count_pocket_contacted_atms
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Count all uniq atoms contacted by the alpha spheres of a given pocket.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pocket *pocket : The pocket
   -----------------------------------------------------------------------------
   ## RETURN:
	int: Number of atoms involed in the pocket
   -----------------------------------------------------------------------------
*/
int count_pocket_contacted_atms(s_pocket *pocket)
{
	if(!pocket || !(pocket->v_lst) || pocket->v_lst->n_vertices <= 0) return -1 ;

	int nb_atoms = 0,
		i = 0 ;

	node_vertice *nvcur = NULL ;
	s_vvertice *vcur = NULL ;

	/* Remember atoms already stored. */
	int atm_ids[pocket->v_lst->n_vertices * 4] ;

	/* Do the search  */
	nvcur = pocket->v_lst->first ;
	while(nvcur) {
		vcur = nvcur->vertice ;
		/*printf("ID in the pocket: %d (%.3f %.3f %.3f\n", vcur->id, vcur->x, vcur->y, vcur->z) ;*/
		for(i = 0 ; i < 4 ; i++) {
			if(in_tab(atm_ids,  nb_atoms, vcur->neigh[i]->id) == 0) {
				atm_ids[nb_atoms] = vcur->neigh[i]->id ;
				nb_atoms ++ ;
			}
		}
		nvcur = nvcur->next ;
	}

	return nb_atoms ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	get_pocket_contacted_atms
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Get pocket vertices under the form of an array of pointer.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pocket *pocket : The pocket
   -----------------------------------------------------------------------------
   ## RETURN:
	s_vvertice**: All pointers to vertices
   -----------------------------------------------------------------------------
*/
s_vvertice** get_pocket_pvertices(s_pocket *pocket)
{
	s_vvertice **pverts = my_calloc(pocket->size, sizeof(s_vvertice*)) ;
	int i = 0 ;
	node_vertice *nvcur = pocket->v_lst->first ;
	while(nvcur) {
		pverts[i] = nvcur->vertice ;
		nvcur = nvcur->next ;
		i++ ;
	}

	return pverts ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	reset_pocket
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Reset pocket descriptors
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pocket *pocket: The pocket.
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/

void reset_pocket(s_pocket *pocket)
{
	pocket->rank = -1 ;

	pocket->score = 0.0 ;
	pocket->ovlp = -2.0 ;
	pocket->ovlp2 = -2.0 ;
	pocket->vol_corresp = -2.0 ;

	pocket->nAlphaPol = 0 ;
	pocket->nAlphaApol = 0 ;

	reset_desc(pocket->pdesc) ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	print_pockets
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Print pockets info in the given buffer
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ FILE *f				 : File to print in
	@ c_lst_pockets *pockets : All pockets
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void print_pockets(FILE *f, c_lst_pockets *pockets) 
{
	node_pocket *cur = pockets->first ;
	if(pockets) {
		fprintf(f, "\n## FPOCKET RESULTS ##\n");
		cur = pockets->first;
		while(cur) {
			print_pocket(f, cur->pocket) ;
			cur = cur->next ;
		}
	}
	else {
		fprintf(f, "\n## NO POCKETS IN THE LIST ##\n");
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	print_pockets_inv
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Print pockets info in the given buffer, starting with the last pocket of the
	chained list.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ FILE *f				: File to print in
	@ c_lst_pockets *pockets: Pockets to print
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void print_pockets_inv(FILE *f, c_lst_pockets *pockets) 
{
	if(pockets) {
		fprintf(f, "\n## POCKET ##\n");
		node_pocket *cur = pockets->last ;
		while(cur) {
			print_pocket(f, cur->pocket) ;
			cur = cur->prev ;
		}
	}
	else {
		fprintf(f, "\n## NO POCKETS IN THE LIST ##\n");
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	void print_pocket(FILE *f, s_pocket *pocket) 
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Print one pocket info in the given buffer
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ FILE *f			: File to print in
	@ s_pocket *pocket  : Pocket to print
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void print_pocket(FILE *f, s_pocket *pocket) 
{
	if(pocket) {
		fprintf(f, "\n## POCKET %d ##\n",pocket->v_lst->first->vertice->resid);
		
		if(pocket->ovlp > -1.0) fprintf(f, "\t Correspondance: %f\n", pocket->ovlp) ;
		if(pocket->vol_corresp > -1.0) fprintf(f, "\t Volume Correspondance: %f\n", pocket->vol_corresp) ;
		fprintf(f, "\t 0  - Pocket Score:                %.4f\n", pocket->score) ;
		fprintf(f, "\t 1  - Number of Voronoi vertices:  %d\n", pocket->pdesc->nb_asph) ;
		fprintf(f, "\t 2  - Mean alpha-sphere radius:    %f\n", pocket->pdesc->mean_asph_ray) ;
		fprintf(f, "\t 3  - Mean alpha-sphere solvent accessibility: %f\n", pocket->pdesc->masph_sacc) ;
		fprintf(f, "\t 4  - Flexibility:                 %f\n", pocket->pdesc->flex) ;
		fprintf(f, "\t 5  - Hydrophobicity Score:        %f\n", pocket->pdesc->hydrophobicity_score) ;
		fprintf(f, "\t 6  - Polarity Score:              %d\n", pocket->pdesc->polarity_score) ;
		fprintf(f, "\t 7  - Volume Score:                %f\n", pocket->pdesc->volume_score) ;
		fprintf(f, "\t 8  - Real volume (approximation): %f\n", pocket->pdesc->volume) ;
		fprintf(f, "\t 9  - Charge Score:                %d\n", pocket->pdesc->charge_score) ;
		fprintf(f, "\t 10 - Local hydrophobic density Score: %f\n", pocket->pdesc->mean_loc_hyd_dens) ;
		fprintf(f, "\t 11 - Number of apolar alpha sphere: %d\n", pocket->nAlphaApol) ;
		fprintf(f, "\t 12 - Amino acid composition:\n") ;
		int i ;
		fprintf(f, "\t 12 - ") ;
		for(i = 0 ; i < M_NB_AA ; i++) fprintf(f, "%3d ", pocket->pdesc->aa_compo[i]) ;
		fprintf(f, "\n") ;
	}
}

