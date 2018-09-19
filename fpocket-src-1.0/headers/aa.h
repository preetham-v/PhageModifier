
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

#ifndef DH_AA
#define DH_AA

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ------------------------------MACROS---------------------------------------*/

#define M_NB_AA 20
#define M_ALA_IDX 0
#define M_CYS_IDX 1
#define M_ASP_IDX 2
#define M_GLU_IDX 3
#define M_PHE_IDX 4
#define M_GLY_IDX 5
#define M_HIS_IDX 6
#define M_ILE_IDX 7
#define M_LYS_IDX 8
#define M_LEU_IDX 9
#define M_MET_IDX 10
#define M_ASN_IDX 11
#define M_PRO_IDX 12
#define M_GLN_IDX 13
#define M_ARG_IDX 14
#define M_SER_IDX 15
#define M_THR_IDX 16
#define M_VAL_IDX 17
#define M_TRP_IDX 18
#define M_TYR_IDX 19

/* --------------------------- PUBLIC STRUCTURES ------------------------------*/

/**
	A structure for the modelisation of an amino acid
*/
typedef struct s_amino_a
{
	char name3[4] ;
	char code ;

	float mw ;
	float volume;
	float hydrophobicity;
	int charge,
	 	polarity, 
		func_grp ;

} s_amino_a ;

/* -----------------------------PROTOTYPES------------------------------------*/

int get_aa_index(const char *name) ; 
char* get_aa_name3(const int index) ;

float get_aa_mw(const char *name) ;
float get_aa_volume_score(const char *name) ;
float get_aa_hydrophobicity_score(const char *name) ;
int get_aa_charge(const char *name)  ;
int get_aa_polarity(const char *name) ;
int get_aa_func_grp(const char *name) ;

float get_volume_score_from_idx(int aa_index) ;
float get_hydrophobicity_score_from_idx(int aa_index) ;
int get_charge_from_idx(int aa_index) ;
int get_polarity_from_idx(int aa_index) ;
int get_func_grp_from_idx(int aa_index) ;

#endif
