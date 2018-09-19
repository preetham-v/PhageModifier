
#include "../headers/writepocket.h"
/**

## ----- GENERAL INFORMATION
##
## FILE 					writepocket.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			02-12-08
##
## ----- SPECIFICATIONS
##
##  Output routine for pockets.
##
## ----- MODIFICATIONS HISTORY
##
##  02-12-08    (v)  Comments UTD
##	01-04-08	(v)  Added template for comments and creation of history
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
	write_single_pdb
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Write atoms and vertices given in argument in the following standard v2.2 
	pdb format.

   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ const char out[]       : Output file name
	@ s_pdb *pdb             : PDB infos
	@ c_lst_pockets *pockets : All pockets
   -----------------------------------------------------------------------------
   ## RETURN:
   -----------------------------------------------------------------------------
*/
void write_pockets_single_pdb(const char out[],  s_pdb *pdb, c_lst_pockets *pockets) 
{
	node_pocket *nextPocket ;
	node_vertice *nextVertice ;
	FILE *f = fopen(out, "w") ;
	if(f) {
		if(pdb) {
			if(pdb->latoms) write_pdb_atoms(f, pdb->latoms, pdb->natoms ) ;
		}

		if(pockets){
			pockets->current = pockets->first ;

			while(pockets->current){
				pockets->current->pocket->v_lst->current = pockets->current->pocket->v_lst->first ;

				while(pockets->current->pocket->v_lst->current){
 					write_pdb_vert(f, pockets->current->pocket->v_lst->current->vertice) ;

					nextVertice = pockets->current->pocket->v_lst->current->next;
					pockets->current->pocket->v_lst->current = nextVertice;
				}

				nextPocket=pockets->current->next;
				pockets->current=nextPocket;
			}
		}

		fclose(f) ;
	}
	else {
		fprintf(stderr, "! The file %s could not be opened!\n", out);
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	write_pdb_atoms
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Print list of atoms as pdb format in given buffer
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ FILE *f      : Buffer to write in.
	@ s_atm *atoms : List of atoms
	@ int natoms   : Number of atoms
   -----------------------------------------------------------------------------
   ## RETURN:
   -----------------------------------------------------------------------------
*/
void write_pdb_atoms(FILE *f, s_atm *atoms, int natoms) 
{
	s_atm *atom = NULL ;
	int i = 0 ;
	for(i = 0 ; i < natoms ; i++) {
		atom = atoms + i ;
		write_pdb_atom_line(f, atom->type, atom->id, atom->name, atom->pdb_aloc, 
							atom->res_name, atom->chain, atom->res_id, 
 							atom->pdb_insert, atom->x, atom->y, atom->z,
 							atom->occupancy, atom->bfactor, atom->symbol, 
 							atom->charge);
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	write_pockets_single_pqr
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Write only pockets (alpha sphere) given in argument in the pqr format.

	!! No atoms writen here, only all pockets in a single pqr file.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ const char out[]       : Output file path
	@ c_lst_pockets *pockets : List of pockets
   -----------------------------------------------------------------------------
   ## RETURN:
   -----------------------------------------------------------------------------
*/
void write_pockets_single_pqr(const char out[], c_lst_pockets *pockets) 
{
	node_pocket *nextPocket ;
	node_vertice *nextVertice ;

	FILE *f = fopen(out, "w") ;
	if(f) {

		if(pockets){
		fprintf(f, "HEADER\n") ;
		fprintf(f, "HEADER This is a pqr format file writen by the programm fpocket.                 \n") ;
		fprintf(f, "HEADER It contains all the pockets vertices found by fpocket.                    \n") ;
			pockets->current = pockets->first ;

			while(pockets->current){
				pockets->current->pocket->v_lst->current = pockets->current->pocket->v_lst->first ;

				while(pockets->current->pocket->v_lst->current){
					write_pqr_vert(f, pockets->current->pocket->v_lst->current->vertice) ;

					nextVertice = pockets->current->pocket->v_lst->current->next;
					pockets->current->pocket->v_lst->current = nextVertice;
				}

				nextPocket=pockets->current->next;
				pockets->current=nextPocket;
			}
		}

		fprintf(f, "TER\nEND\n") ;
		fclose(f) ;
	}
	else {
		fprintf(stderr, "! The file %s could not be opened!\n", out);
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	write_each_pocket
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Write each pocket in a single pqr (vertices) and pdb (atoms) file format.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ const char out[]       : Output file path
	@ c_lst_pockets *pockets : List of pockets
   -----------------------------------------------------------------------------
   ## RETURN:
   -----------------------------------------------------------------------------
*/
void write_each_pocket(const char out_path[], c_lst_pockets *pockets)
{
	int out_len = strlen(out_path) ;
	char out[out_len+20] ;
	out[0] = '\0' ;

	node_pocket *pcur ;
	
	int i = 0 ;
	if(pockets){
		pcur = pockets->first ;

		while(pcur){
			sprintf(out, "%s/pocket%d_vert.pqr", out_path, i) ;
			write_pocket_pqr(out, pcur->pocket) ;

			sprintf(out, "%s/pocket%d_atm.pdb", out_path, i) ;
			write_pocket_pdb(out, pcur->pocket) ;

			pcur = pcur->next ;
			i++ ;
		}
	}
	else {
		fprintf(stderr, "! The file %s could not be opened!\n", out);
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	void write_pocket_pqr
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Write vertices of the pocket given in argument in the pqr format.

   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ const char out[]  : Output file path
	@ s_pocket *pocket : The pocket to write
   -----------------------------------------------------------------------------
   ## RETURN:
   -----------------------------------------------------------------------------
*/
void write_pocket_pqr(const char out[], s_pocket *pocket) 
{
	node_vertice *vcur = NULL ;

	FILE *f = fopen(out, "w") ;
	if(f && pocket) {
		fprintf(f, "HEADER\n") ;
		fprintf(f, "HEADER This is a pqr format file writen by the programm fpocket.                 \n") ;
		fprintf(f, "HEADER It represent the voronoi vertices of a single pocket found by the         \n") ;
		fprintf(f, "HEADER algorithm.                                                                \n") ;
		fprintf(f, "HEADER                                                                           \n") ;
		fprintf(f, "HEADER Information about the pocket %5d:\n", pocket->v_lst->first->vertice->resid) ;
		fprintf(f, "HEADER 0  - Pocket Score                      : %.4f\n", pocket->score) ;
		fprintf(f, "HEADER 1  - Number of V. Vertices             : %5d\n", pocket->pdesc->nb_asph) ;
		fprintf(f, "HEADER 2  - Mean alpha-sphere radius          : %.4f\n", pocket->pdesc->mean_asph_ray) ;
		fprintf(f, "HEADER 3  - Mean alpha-sphere SA              : %.4f\n", pocket->pdesc->masph_sacc) ;
		fprintf(f, "HEADER 4  - Mean B-factor                     : %.4f\n", pocket->pdesc->flex) ;
		fprintf(f, "HEADER 5  - Hydrophobicity Score              : %.4f\n", pocket->pdesc->hydrophobicity_score) ;
		fprintf(f, "HEADER 6  - Polarity Score                    : %5d\n", pocket->pdesc->polarity_score) ;
		fprintf(f, "HEADER 7  - Volume Score                      : %.4f\n", pocket->pdesc->volume_score) ;
		fprintf(f, "HEADER 8  - Real volume (approximation)       : %.4f\n", pocket->pdesc->volume) ;
		fprintf(f, "HEADER 9  - Charge Score                      : %5d\n", pocket->pdesc->charge_score) ;
		fprintf(f, "HEADER 10 - Local hydrophobic density Score   : %.4f\n", pocket->pdesc->mean_loc_hyd_dens) ;
		fprintf(f, "HEADER 11 - Number of apolar alpha sphere     : %5d\n", pocket->nAlphaApol) ;
		fprintf(f, "HEADER 12 - Proportion of apolar alpha sphere : %.4f\n", pocket->pdesc->apolar_asphere_prop) ;

		vcur = pocket->v_lst->first ;

		while(vcur){
			write_pqr_vert(f, vcur->vertice) ;

			vcur = vcur->next ;
		}

		fprintf(f, "TER\nEND\n") ;
		fclose(f) ;
	}
	else {
		if(!f) fprintf(stderr, "! The file %s could not be opened!\n", out);
		else fprintf(stderr, "! Invalid pocket to write in write_pocket_pqr !\n");
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	write_pocket_pdb
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Write atoms contacted by vertices of the pocket given in argument in 
	the pdb format.

   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ const char out[]  : Output file path
	@  s_pocket *pocket : The pocket to write
   -----------------------------------------------------------------------------
   ## RETURN:
   -----------------------------------------------------------------------------
*/
void write_pocket_pdb(const char out[], s_pocket *pocket) 
{
	node_vertice *vcur = NULL ;
	int i = 0 ;
	int cur_size = 0,
		cur_allocated = 10 ;

	s_atm **atms = (s_atm**)my_malloc(sizeof(s_atm*)*10) ;
	s_atm *atom = NULL ;

	FILE *f = fopen(out, "w") ;
	if(f && pocket) {

		fprintf(f, "HEADER\n") ;
		fprintf(f, "HEADER This is a pdb format file writen by the programm fpocket.                 \n") ;
		fprintf(f, "HEADER It represents the atoms contacted by the voronoi vertices of the pocket.  \n") ;
		fprintf(f, "HEADER                                                                           \n") ;
		fprintf(f, "HEADER Information about the pocket %5d:\n", pocket->v_lst->first->vertice->resid) ;
		fprintf(f, "HEADER 0  - Pocket Score                      : %.4f\n", pocket->score) ;
		fprintf(f, "HEADER 1  - Number of V. Vertices             : %5d\n", pocket->pdesc->nb_asph) ;
		fprintf(f, "HEADER 2  - Mean alpha-sphere radius          : %.4f\n", pocket->pdesc->mean_asph_ray) ;
		fprintf(f, "HEADER 3  - Mean alpha-sphere SA              : %.4f\n", pocket->pdesc->masph_sacc) ;
		fprintf(f, "HEADER 4  - Mean B-factor                     : %.4f\n", pocket->pdesc->flex) ;
		fprintf(f, "HEADER 5  - Hydrophobicity Score              : %.4f\n", pocket->pdesc->hydrophobicity_score) ;
		fprintf(f, "HEADER 6  - Polarity Score                    : %5d\n", pocket->pdesc->polarity_score) ;
		fprintf(f, "HEADER 7  - Volume Score                      : %.4f\n", pocket->pdesc->volume_score) ;
		fprintf(f, "HEADER 8  - Real volume (approximation)       : %.4f\n", pocket->pdesc->volume) ;
		fprintf(f, "HEADER 9  - Charge Score                      : %5d\n", pocket->pdesc->charge_score) ;
		fprintf(f, "HEADER 10 - Local hydrophobic density Score   : %.4f\n", pocket->pdesc->mean_loc_hyd_dens) ;
		fprintf(f, "HEADER 11 - Number of apolar alpha sphere     : %5d\n", pocket->nAlphaApol) ;
		fprintf(f, "HEADER 12 - Proportion of apolar alpha sphere : %.4f\n", pocket->pdesc->apolar_asphere_prop) ;

	/* First get the list of atoms */
		vcur = pocket->v_lst->first ;

		while(vcur){
			for(i = 0 ; i < 4 ; i++) {
				if(!is_in_lst_atm(atms, cur_size, vcur->vertice->neigh[i]->id)) {
					if(cur_size >= cur_allocated-1) {
						cur_allocated *= 2 ;
						atms = (s_atm**)my_realloc(atms, sizeof(s_atm)*cur_allocated) ;
					}
					atms[cur_size] = vcur->vertice->neigh[i] ;
					cur_size ++ ;
				}

			}
			vcur = vcur->next ;
		}

	/* Then write atoms... */

		for(i = 0 ; i < cur_size ; i++) {
			atom = atms[i] ;

			write_pdb_atom_line(f, atom->type, atom->id, atom->name, atom->pdb_aloc, 
									atom->res_name, atom->chain, atom->res_id, 
 									atom->pdb_insert, atom->x, atom->y, atom->z,
 									atom->occupancy, atom->bfactor, atom->symbol, 
 									atom->charge);
		}
	
		fprintf(f, "TER\nEND\n") ;
		fclose(f) ;
	}
	else {
		if(!f) fprintf(stderr, "! The file %s could not be opened!\n", out);
		else fprintf(stderr, "! Invalid pocket to write in write_pocket_pqr !\n");
	}

	my_free(atms) ;
}
