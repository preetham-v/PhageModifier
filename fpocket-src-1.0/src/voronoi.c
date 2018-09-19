
#include "../headers/voronoi.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					voronoi.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			02-12-08
##
## ----- SPECIFICATIONS
##
##  Functions dealing with input/output of qhull: we will send to qhull a set
##  of points (basically all atoms of the system), and qhull will send back all
##  voronoi vertices positions. Everything will then be parsed properly for
##  future use.
##
## ----- MODIFICATIONS HISTORY
##
##	02-12-08	(v)  Comments UTD
##	01-04-08	(v)  Added template for comments and creation of history
##	21-02-08	(p)	 Adding support for proteins with hydrogens
##	01-01-08	(vp) Created (random date...)
##
## ----- TODO or SUGGESTIONS
##
##  (v) Improve the reading of vertices. Better: include qhull algorithm in our
##		source code.
##  (v) Handle errors when reading/filling vertices..........

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

static void fill_vvertices(s_lst_vvertice *lvvert, const char fpath[], s_atm *atoms, int natoms,
						   int min_apol_neigh, float asph_min_size, float asph_max_size) ;

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	s_lst_vvertice
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Calculate voronoi vertices using an ensemble of atoms, and then load resulting
	vertices into a s_lst_vvertice structure. The function call an external
	programm qvoronoi, part of qhull programme which can be download at:
		http://www.qhull.org/download/
	or installed with apt-get install qhull-bin
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pdb *pdb          : PDB informations
	@ int min_apol_neigh  : Number of apolar neighbor of a vertice to be
							considered as apolar
	@ float asph_min_size : Minimum size of voronoi vertices to retain
	@ float asph_max_size : Maximum size of voronoi vertices to retain
   -----------------------------------------------------------------------------
   ## RETURN:
	s_lst_vvertice * :The structure containing the list of vertices.
   -----------------------------------------------------------------------------
*/
s_lst_vvertice* load_vvertices(s_pdb *pdb, int min_apol_neigh, float asph_min_size, float asph_max_size)
{
	int i, nb_h=0;
	s_atm *ca = NULL ;
	s_lst_vvertice *lvvert = NULL ;
    FILE *ftmp=fopen("/tmp/fpocket_qvor.dat","w");
	FILE *fvoro = fopen("/tmp/qvoro_tmp.dat", "w+");
/* 	lvvert->h_tr=(int *)my_malloc(sizeof(int));*/



	if(fvoro != NULL) {
		lvvert = (s_lst_vvertice *)my_malloc(sizeof(s_lst_vvertice)) ;
		lvvert->h_tr=NULL;
		/* Loop a first time to get out how many heavy atoms are in the file */
		for(i = 0; i <  pdb->natoms ; i++){
			ca = (pdb->latoms)+i ;
			if(strcmp(ca->symbol,"H")) {
				lvvert->h_tr=(int *)my_realloc(lvvert->h_tr,sizeof(int)*(i-nb_h+1)) ;
				lvvert->h_tr[i-nb_h]=i ;
			}
			else nb_h++;
		}
		lvvert->n_h_tr=i-nb_h;

		/* Write the header for qvoronoi */
		fprintf(fvoro,"3 rbox D3\n%d\n", lvvert->n_h_tr);
		/* Loop a second time for the qvoronoi input coordinates */
		for(i = 0; i <  pdb->natoms ; i++){
			ca = (pdb->latoms)+i ;
			if(strcmp(ca->symbol,"H")) {
			/* Only if this is a heavy atom export it for voronoi tesselation,
			 * else discard it */
				fprintf(fvoro,"%f %f %f \n", ca->x, ca->y, ca->z);
			}
		}

		fflush(fvoro) ;
		rewind(fvoro);
		
		//int status = system("qvoronoi p i Pp Fn < voro_tmp.dat > voro.tmp") ;
        run_qvoronoi(fvoro,ftmp);
        int status=M_VORONOI_SUCCESS;
		
		fclose(fvoro) ;
        fclose(ftmp);
		
        remove("/tmp/qvoro_tmp.dat");
		if(status == M_VORONOI_SUCCESS) {
			fill_vvertices(lvvert, "/tmp/fpocket_qvor.dat", pdb->latoms, pdb->natoms,
							min_apol_neigh, asph_min_size, asph_max_size);
		}
		else {
			my_free(lvvert);
			lvvert = NULL ;
			fprintf(stderr, "! Voronoi command failed with status %d...\n", status) ;
		}
        remove("/tmp/fpocket_qvor.dat");
	}
	else {
		fprintf(stderr, "! File for Voronoi vertices calculation couldn't be opened...\n") ;
	}
        
	return lvvert ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	fill_vertices
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_lst_vvertice *lvvert : The structure to fill
	@ const char fpath[]     : File containing vertices
	@ s_atm *atoms           : List of atoms
	@ int natoms             : Number of atoms
	@ int min_apol_neigh  : Number of apolar neighbor of a vertice to be
							considered as apolar
	@ float asph_min_size : Minimum size of voronoi vertices to retain
	@ float asph_max_size : Maximum size of voronoi vertices to retain
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Fill structure given in argument (must have been allocated) using a file
	containing vertice coordinates and neighbours using p i options of qhull.
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
static void fill_vvertices(s_lst_vvertice *lvvert, const char fpath[], s_atm *atoms, int natoms,
						   int min_apol_neigh, float asph_min_size, float asph_max_size)
{
	FILE *f = NULL ;	/* File handler for vertices coordinates */
	FILE *fNb = NULL ;	/* File handler for vertices atomic neighbours */
	FILE *fvNb = NULL;	/* File handler for vertices vertice neighbours */

	s_vvertice *v = NULL ;

	float tmpRay;	/* Temporary Ray of voronoi vertice (ray of alpha sphere) */
	float xyz[3] = {0,0,0};

	int i, j,nchar_max = 255,
		vInMem = 0,				/* Saved vertices counter */
		curLineNb = 0,			/* Current line number */
		trash = 0,
		tmpApolar=0,
		curVnbIdx[4],			/* Current vertice neighbours */
		curNbIdx[4];			/* Current atomic neighbours */

	char cline[nchar_max],
		 nbline[nchar_max],
		 vNbline[nchar_max],
		 *s_nvert = (char *) my_malloc(sizeof(char)*nchar_max) ;

/* Once we have the number of lines, lets allocate memory and get the lines */
       
	f = fopen(fpath,"r") ;
	fNb = fopen(fpath,"r") ;
	fvNb = fopen(fpath,"r") ;

	char *status = NULL ;

	status = fgets(cline, nchar_max, f) ;		/* Skip fir=st line */
	status = fgets(cline, nchar_max, f) ;		/* Load 2nd line containing nbr of coors. */
	status = fgets(nbline, nchar_max, fNb) ;	/* Skip first line */
	status = fgets(nbline, nchar_max, fNb) ;	/* Load 2nd line containing nbr of coors. */
	status = fgets(vNbline, nchar_max, fvNb) ;	/* Skip first line */
	status = fgets(vNbline, nchar_max, fvNb) ;	/* Load 2nd line containing nbr of coors. */

 	sscanf(cline,"%d",&(lvvert->nvert)) ;
	lvvert->qhullSize = lvvert->nvert ;
	lvvert->tr = (int *) my_malloc(lvvert->nvert*sizeof(int));
	for(i = 0 ; i < lvvert->nvert ; i++) lvvert->tr[i] = -1;

 	lvvert->vertices = (s_vvertice *) my_calloc(lvvert->nvert, sizeof(s_vvertice)) ;
	lvvert->pvertices= (s_vvertice **) my_calloc(lvvert->nvert, sizeof(s_vvertice*)) ;
	
	/* Get the string of number of vertices to read, to look up the neighbour
	 * list from qhull */
	sprintf(s_nvert,"%d",lvvert->nvert);
	strcat(s_nvert,"\n");

	/* Advance cursor to neighbour list */
	while(fgets(nbline, nchar_max, fNb) != NULL && strcmp(s_nvert, nbline) != 0) ;
	/* Advance cursor to the vertice neighbour list */
	while(fgets(vNbline,nchar_max,fvNb) != NULL && curLineNb++ < lvvert->nvert*2+1) ;

	i = 0 ;
	while(fgets(cline, nchar_max, f) != NULL) {
	/* Read vertice positions */
 		if(fgets(nbline, nchar_max,fNb)!=NULL){
		/* Read neighbours */
			if(fgets(vNbline, nchar_max,fvNb)!=NULL){
			/* Read vertice neighbour vertices */
				if(strcmp("\n", cline) != 0 && strcmp("\n", nbline) != 0
				   && strcmp("\n", vNbline) != 0) {

					sscanf(cline,"%f %f %f",&xyz[0], &xyz[1], &xyz[2]);
					sscanf(nbline,"%d %d %d %d",&curNbIdx[0],&curNbIdx[1],
												&curNbIdx[2],&curNbIdx[3]);
 					sscanf(vNbline,"%d %d %d %d %d", &trash, &curVnbIdx[0],
													 &curVnbIdx[1],
													 &curVnbIdx[2], &curVnbIdx[3]);
				/* Test voro. vert. for alpha sphere cond. and returns ray if
				 * cond. are ok, -1 else */
					tmpRay = testVvertice(xyz, curNbIdx, atoms, asph_min_size,
										  asph_max_size,lvvert);
					if(tmpRay > 0){
						v = (lvvert->vertices + vInMem) ;
						v->x = xyz[0]; v->y = xyz[1]; v->z = xyz[2];
						v->ray = tmpRay;
						v->sort_x = -1 ;
						v->seen = 0 ;
						
						tmpApolar=0;

						for(j = 0 ; j < 4 ; j++) {
							v->neigh[j] = &(atoms[lvvert->h_tr[curNbIdx[j]]]);

							if(atoms[lvvert->h_tr[curNbIdx[j]]].electroneg<2.8) tmpApolar++ ;
							if(curVnbIdx[j]>0) v->vneigh[j] = curVnbIdx[j];
						}

						v->apol_neighbours = 0 ;
						lvvert->tr[i] = vInMem ;

						lvvert->pvertices[vInMem] = v ;

						vInMem++ ;		/* Vertices actually read */
						v->id = natoms+i+1-vInMem ;

						if(tmpApolar >= min_apol_neigh) v->type = M_APOLAR_AS;
						else v->type = M_POLAR_AS;

						v->qhullId = i;		/* Set index in the qhull file */
						v->resid = -1;		/* Initialize internal index */
						set_barycenter(v) ;	/* Set barycentre */
					}
					i++ ;
				}
			}
 		}
	}

	lvvert->nvert=vInMem ;
	fclose(f) ;
	fclose(fNb) ;
	fclose(fvNb);
        remove("/tmp/fpocket_qvor.dat");
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	set_barycenter
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Set barycenter of a vertice using it's 4 contacting atoms.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ s_vvertice *v: The vertice
   -----------------------------------------------------------------------------
   ## RETURN: void
   -----------------------------------------------------------------------------
*/
void set_barycenter(s_vvertice *v)
{
	int i ;
	float xsum = 0.0,
		  ysum = 0.0,
		  zsum = 0.0 ;

	for(i = 0 ; i < 4 ; i++) {
		xsum += v->neigh[i]->x ;
		ysum += v->neigh[i]->y ;
		zsum += v->neigh[i]->z ;
	}

	v->bary[0] = xsum*0.25 ;
	v->bary[1] = ysum*0.25 ;
	v->bary[2] = zsum*0.25 ;

}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	testVvertice
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Test if alpha sphere conditions are fulfilled for current vertice
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ float xyz[3]        : Coordinates of current vertice
	@ int curNbIdx[4]     : Indexes of atomic neighbours of the current vertice
	@ s_atm *atoms        : List of all atoms
	@ float min_asph_size : Minimum size of alpha spheres.
	@ float max_asph_size : Maximum size of alpha spheres.
   -----------------------------------------------------------------------------
   ## RETURN:
	float : -1 if conditions are not fulfilled, else the alpha sphere radius
		    is returned.
   -----------------------------------------------------------------------------
*/
float testVvertice(float xyz[3], int curNbIdx[4], s_atm *atoms,
				   float min_asph_size, float max_asph_size,
				   s_lst_vvertice *lvvert)
{
	float x = xyz[0],
		  y = xyz[1],
		  z = xyz[2] ;

	s_atm *cura = &(atoms[lvvert->h_tr[curNbIdx[0]]]) ;

	float distVatom1 = dist(x, y, z, cura->x, cura->y, cura->z) ;
	float distVatom2,
		  distVatom3,
		  distVatom4;

	if(min_asph_size <= distVatom1  && distVatom1 <= max_asph_size){
		cura = &(atoms[lvvert->h_tr[curNbIdx[1]]]) ;
		distVatom2 = dist(x, y, z, cura->x, cura->y, cura->z);

		cura = &(atoms[lvvert->h_tr[curNbIdx[2]]]) ;
 		distVatom3 = dist(x, y, z, cura->x, cura->y, cura->z);

		cura = &(atoms[lvvert->h_tr[curNbIdx[3]]]) ;
  		distVatom4=dist(x, y, z, cura->x, cura->y, cura->z);

		/* Test if all 4 neighbours are on the alpha sphere surface
		 * (approximate test) */
		if(fabs(distVatom1-distVatom2) < M_PREC_TOLERANCE &&
		   fabs(distVatom1-distVatom3)  < M_PREC_TOLERANCE &&
		   fabs(distVatom1-distVatom4) < M_PREC_TOLERANCE){

			return distVatom1;
		}

	}
	return -1.0;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	print_vvertices
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Print function.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ FILE *f                : Buffer to print in
	@ s_lst_vvertice *lvvert : Vertices to print
   -----------------------------------------------------------------------------
   ## RETURN: void
   -----------------------------------------------------------------------------
*/
void print_vvertices(FILE *f, s_lst_vvertice *lvvert)
{
	if(lvvert) {
		if(lvvert->vertices) {
			int i ;
			for(i = 0 ; i < lvvert->nvert ; i++) {
				s_vvertice *v = &(lvvert->vertices[i]) ;
				if( v->neigh[0] &&  v->neigh[1] && v->neigh[2] &&  v->neigh[3]) {
					fprintf(f, "====== Vertice %d: =====\n", i);
					fprintf(f, "- x = %f, y = %f, z = %f\n", v->x, v->y, v->z);
					fprintf(f, "- ix = %d\n",v->sort_x);

					float d1 = dist(v->x, v->y, v->z, v->neigh[0]->x, v->neigh[0]->y, v->neigh[0]->z) ;
					float d2 = dist(v->x, v->y, v->z, v->neigh[1]->x, v->neigh[1]->y, v->neigh[1]->z) ;
					float d3 = dist(v->x, v->y, v->z, v->neigh[2]->x, v->neigh[2]->y, v->neigh[2]->z) ;
					float d4 = dist(v->x, v->y, v->z, v->neigh[3]->x, v->neigh[3]->y, v->neigh[3]->z) ;

					fprintf(f, "- Neighbour: \n1 - %f (%f %f %f: %d) \n2 - %f (%f %f %f: %d)\n3 - %f (%f %f %f: %d)\n4 - %f (%f %f %f: %d)\n", d1, v->neigh[0]->x, v->neigh[0]->y, v->neigh[0]->z, v->neigh[0]->id, d2, v->neigh[1]->x, v->neigh[1]->y, v->neigh[1]->z, v->neigh[1]->id, d3, v->neigh[2]->x,  v->neigh[2]->y ,  v->neigh[2]->z, v->neigh[2]->id, d4, v->neigh[3]->x, v->neigh[3]->y, v->neigh[3]->z, v->neigh[3]->id);
				}
			}
		}
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	get_verts_volume_ptr
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Get an monte carlo approximation of the volume occupied by the alpha spheres
	given in argument (list of pointers)
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_vvertice **verts: List of pointer to alpha spheres
	@ int nvert: Number of spheres
	@ int niter: Number of monte carlo iteration to perform
   -----------------------------------------------------------------------------
   ## RETURN:
	float: volume.
   -----------------------------------------------------------------------------
*/
float get_verts_volume_ptr(s_vvertice **verts, int nvert, int niter)
{
	int i = 0, j = 0,
		nb_in = 0;

	float xmin = 0.0, xmax = 0.0,
		  ymin = 0.0, ymax = 0.0,
		  zmin = 0.0, zmax = 0.0,
		  xtmp = 0.0, ytmp = 0.0, ztmp = 0.0,
		  xr   = 0.0, yr   = 0.0, zr   = 0.0,
		  vbox = 0.0 ;

	s_vvertice *vcur = NULL ;

	/* First, search extrems coordinates to get a contour box of the molecule */
	for(i = 0 ; i < nvert ; i++) {
		vcur = verts[i] ;

		if(i == 0) {
			xmin = vcur->x - vcur->ray ; xmax = vcur->x + vcur->ray ;
			ymin = vcur->y - vcur->ray ; ymax = vcur->y + vcur->ray ;
			zmin = vcur->z - vcur->ray ; zmax = vcur->z + vcur->ray ;
		}
		else {
		/* Update the minimum and maximum extreme point */
			if(xmin > (xtmp = vcur->x - vcur->ray)) xmin = xtmp ;
			else if(xmax < (xtmp = vcur->x + vcur->ray)) xmax = xtmp ;

			if(ymin > (ytmp = vcur->y - vcur->ray)) ymin = ytmp ;
			else if(ymax < (ytmp = vcur->y + vcur->ray)) ymax = ytmp ;

			if(zmin > (ztmp = vcur->z - vcur->ray)) zmin = ztmp ;
			else if(zmax < (ztmp = vcur->z + vcur->ray)) zmax = ztmp ;
		}
	}

	/* Next calculate the contour box volume */
	vbox = (xmax - xmin)*(ymax - ymin)*(zmax - zmin) ;

	/* Then apply monte carlo approximation of the volume.	 */
	for(i = 0 ; i < niter ; i++) {
		xr = rand_uniform(xmin, xmax) ;
		yr = rand_uniform(ymin, ymax) ;
		zr = rand_uniform(zmin, zmax) ;

		for(j = 0 ; j < nvert ; j++) {
			vcur = verts[j] ;
			xtmp = vcur->x - xr ;
			ytmp = vcur->y - yr ;
			ztmp = vcur->z - zr ;

		/* Compare r^2 and dist(center, random_point)^2 */
			if((vcur->ray*vcur->ray) > (xtmp*xtmp + ytmp*ytmp + ztmp*ztmp)) {
			/* The point is inside one of the vertice!! */
				nb_in ++ ; break ;
			}
		}
	}

	/* Ok lets just return the volume Vpok = Nb_in/Niter*Vbox */
	return ((float)nb_in)/((float)niter)*vbox;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	free_vert_lst
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Free memory
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ s_lst_vvertice *lvvert : Data to free
   -----------------------------------------------------------------------------
   ## RETURN: void
   -----------------------------------------------------------------------------
*/
void free_vert_lst(s_lst_vvertice *lvvert)
{
	if(lvvert) {
		if(lvvert->vertices) {
			my_free(lvvert->vertices) ;
			lvvert->vertices = NULL ;
		}
		if(lvvert->pvertices) {
			my_free(lvvert->pvertices) ;
			lvvert->pvertices = NULL ;
		}
		if(lvvert->tr) {
			my_free(lvvert->tr) ;
			lvvert->tr = NULL ;
		}
		if(lvvert->h_tr) {
			my_free(lvvert->h_tr) ;
			lvvert->h_tr = NULL ;
		}
		my_free(lvvert) ;
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	is_in_lst_vert
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Says wether a vertice of id v_id is in a list of vertices or not
   -----------------------------------------------------------------------------
   ## PARAMETRES:
   -----------------------------------------------------------------------------
   ## RETURN:
	1 if the vertice is in the tab, 0 if not
   -----------------------------------------------------------------------------
*/
int is_in_lst_vert(s_vvertice **lst_vert, int nb_vert, int v_id)
{
	int i ;
	for(i = 0 ; i < nb_vert ; i++) {
		if(v_id == lst_vert[i]->id) return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	is_in_lst_vert
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Says wether a vertice of id v_id is in a list of vertices or not
   -----------------------------------------------------------------------------
   ## PARAMETRES:
   -----------------------------------------------------------------------------
   ## RETURN:
	1 if the vertice is in the tab, 0 if not
   -----------------------------------------------------------------------------
*/
int is_in_lst_vert_p(s_vvertice **lst_vert, int nb_vert, s_vvertice *vert)
{
	int i ;
	for(i = 0 ; i < nb_vert ; i++) {
		if(vert == lst_vert[i]) return 1 ;
	}

	return 0 ;
}
/** -----------------------------------------------------------------------------
	-----------------------------------------------------------------------------
	-----------------------------------------------------------------------------

	OUTPUT FUNCTIONS

	-----------------------------------------------------------------------------
	-----------------------------------------------------------------------------
	-----------------------------------------------------------------------------
*/

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	void write_pdb_vertice(FILE *f, s_vvertice *v)
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Write a voronoi vertice in pdb format.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ FILE *f: file to write in
	@ s_vvertice *v: The vertice
   -----------------------------------------------------------------------------
   ## RETURN:
   -----------------------------------------------------------------------------
*/
void write_pdb_vert(FILE *f, s_vvertice *v)
{
	if(v->type==M_APOLAR_AS) write_pdb_atom_line(f, "HETATM", v->id, "APOL",
												 ' ', "STP", "C", v->resid, ' ',
												 v->x, v->y, v->z, 0.0, 0.0,
												 "Ve", -1);

	else write_pdb_atom_line(f, "HETATM", v->id, " POL", ' ', "STP", "C",
								v->resid, ' ',v->x, v->y, v->z,0.0, 0.0,
								"Ve", -1);
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	write_pqr_vertice
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Write a voronoi vertice in pqr format.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ FILE *f       : file to write in
	@ s_vvertice *v : The vertice
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void write_pqr_vert(FILE *f, s_vvertice *v)
{
	if(v->type==M_APOLAR_AS) write_pqr_atom_line(f, "ATOM", v->id, "APOL", ' ',
												 "STP", " ", v->resid, ' ',
												  v->x, v->y, v->z, 0.0, v->ray);

	else write_pqr_atom_line(f, "ATOM", v->id, " POL", ' ', "STP", " ",
							 v->resid, ' ',v->x, v->y, v->z,0.0, v->ray);
}

