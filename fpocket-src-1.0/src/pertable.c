
#include "../headers/pertable.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					pertable.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			28-11-08
##
## ----- SPECIFICATIONS
##
## This file defines the periodic element table. It's strongly based on the
## VMD source code.
##
## ----- MODIFICATIONS HISTORY
##
##  17-03-09    (v)  Added function testing if a string is a valid element symbol
##	28-11-08	(v)  Comments UTD
##	01-04-08	(v)  Added template for comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## ----- TODO or SUGGESTIONS
##
##	(v) Merge with atom.c
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

static const int ST_nelem = 112 ;

static const char *ST_pte_symbol[] = { 
	"X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
	"Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
	"Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", 
	"As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
	"Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
	"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
	"Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
	"Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
	"Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
	"Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
	"Ds", "Rg"
} ;

static const float ST_pte_electronegativity[] = {
	 0.0,  2.1, 0.98,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0, -1.0,
	 0.9,  1.2,  1.5,  1.8,  2.1,  2.5,  3.0, -1.0,  0.8,  1.0,  1.3,
	 1.5,  1.6,  1.6,  1.5,  1.8,  1.8,  1.9,  1.9,  1.6,  1.8,  2.0,
	 2.2,  2.4,  2.9, -1.0,  0.8,  1.0,  1.2,  1.3,  1.6,  2.0,  1.9,
	 2.2,  2.2,  2.3,  1.9,  1.7,  1.7,  1.8,  2.0,  2.1,  2.6,  2.6,
	 0.8,  0.9, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
	-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  1.3,  1.5,  1.7,  1.9,  2.2,
	 2.2,  2.2,  2.4,  1.9,  1.8,  1.8,  1.9,  2.0,  2.2, -1.0,  0.7,
	 0.9,  1.1,  1.3,  1.5,  1.7,  1.3,  1.3,  1.3,  1.3,  1.3,  1.3,
	 1.3,  1.3,  1.3,  1.3, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
	-1.0, -1.0
} ;

static const float ST_pte_mass[] = { 
	/* X  */ 0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811,  
	/* C  */ 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797, 
	/* Na */ 22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
	/* S  */ 32.065, 35.453, 39.948, 39.0983, 40.078, 44.955910,
	/* Ti */ 47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332,
	/* Ni */ 58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160, 
	/* Se */ 78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585, 
	/* Zr */ 91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550,
	/* Pd */ 106.42, 107.8682, 112.411, 114.818, 118.710, 121.760, 
	/* Te */ 127.60, 126.90447, 131.293, 132.90545, 137.327, 
	/* La */ 138.9055, 140.116, 140.90765, 144.24, 145.0, 150.36,
	/* Eu */ 151.964, 157.25, 158.92534, 162.500, 164.93032, 
	/* Er */ 167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,
	/* W  */ 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655, 
	/* Hg */ 200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0, 
	/* Fr */ 223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,
	/* Np */ 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
	/* Md */ 258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,
	/* Mt */ 268.0, 271.0, 272.0
};

/*
	 A. Bondi, J. Phys. Chem., 68, 441 - 452, 1964, 
	.Phys.Chem., 100, 7384 - 7391, 1996.
 */
static const float ST_pte_rvdw[] = {
	/* X  */ 1.5, 1.2, 1.4, 1.82, 2.0, 2.0,  
	/* C  */ 1.7, 1.55, 1.52, 1.47, 1.54, 
	/* Na */ 2.27, 1.73, 2.0, 2.1, 1.8,
	/* S  */ 1.8, 1.75, 1.88, 2.75, 2.0, 2.0,
	/* Ti */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
	/* Ni */ 1.63, 1.4, 1.39, 1.07, 2.0, 1.85,
	/* Se */ 1.9, 1.85, 2.02, 2.0, 2.0, 2.0, 
	/* Zr */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
	/* Pd */ 1.63, 1.72, 1.58, 1.93, 2.17, 2.0, 
	/* Te */ 2.06, 1.98, 2.16, 2.0, 2.0,
	/* La */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
	/* Eu */ 2.0, 2.0, 2.0, 2.0, 2.0,
	/* Er */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
	/* W  */ 2.0, 2.0, 2.0, 2.0, 1.72, 1.66,
	/* Hg */ 1.55, 1.96, 2.02, 2.0, 2.0, 2.0, 2.0,
	/* Fr */ 2.0, 2.0, 2.0, 2.0, 2.0, 1.86,
	/* Np */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
	/* Md */ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
	/* Mt */ 2.0, 2.0, 2.0
};

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	pte_get_mass
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Returns the mass for a given element
   -----------------------------------------------------------------------------
   ## PARAMETERS:	
	@ const char *symbol: The symbol of the element in the periodic table
   -----------------------------------------------------------------------------
   ## RETURN:
	float: mass corresponding to symbol
   -----------------------------------------------------------------------------
*/
float pte_get_mass(const char *symbol)
{	
	char atom[3] ;
	if (symbol != NULL) {
		atom[0] = (char) toupper((int) symbol[0]);
		atom[1] = (char) tolower((int) symbol[1]);	
		atom[2] = '\0' ;
	
		int i ;
		for (i = 0; i < ST_nelem ; i++) {
			if ( (ST_pte_symbol[i][0] == atom[0]) && (ST_pte_symbol[i][1] == atom[1]) ) {
				
				return ST_pte_mass[i] ;
			}
		}
	}


	return -1 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	pte_get_vdw_ray
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Returns the van der walls radius for a given element
   -----------------------------------------------------------------------------
   ## PARAMETERS:	
	@ const char *symbol: The symbol of the element in the periodic table
   -----------------------------------------------------------------------------
   ## RETURN:
	float: vdw radius corresponding to symbol
   -----------------------------------------------------------------------------
*/
float pte_get_vdw_ray(const char *symbol)
{
	char atom[3] ;

	if (symbol != NULL) {
		atom[0] = (char) toupper((int) symbol[0]);
		atom[1] = (char) tolower((int) symbol[1]);
		atom[2] = '\0' ;
	
		int i ;
		for (i = 0; i < ST_nelem ; i++) {
			if ( (ST_pte_symbol[i][0] == atom[0]) && (ST_pte_symbol[i][1] == atom[1]) ) {
				return ST_pte_rvdw[i] ;
			}
		}
	}

	return -1 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	pte_get_enegativity
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Returns the electronegativity (Pauling) value for a given element
   -----------------------------------------------------------------------------
   ## PARAMETERS:	
	@ const char *symbol: The symbol of the element in the periodic table
   -----------------------------------------------------------------------------
   ## RETURN:
	float: electrobegativity of Pauling corresponding to symbol
   -----------------------------------------------------------------------------
*/
float pte_get_enegativity(const char *symbol)
{
	char atom[3] = "" ;

	if (symbol != NULL) {
		atom[0] = (char) toupper((int) symbol[0]);
		atom[1] = (char) tolower((int) symbol[1]);
		atom[2] = '\0' ;
	
		int i ;
		for (i = 0; i < ST_nelem ; i++) {
			if ( (ST_pte_symbol[i][0] == atom[0]) && (ST_pte_symbol[i][1] == atom[1]) ) {
				return ST_pte_electronegativity[i] ;
			}
		}
	}

	return -1 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	is_valid_element
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Check if a given string corresponds to an atom element.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ const char *str : The string to test
	@ int tcase       : If = 1, dont take into account the case.
   -----------------------------------------------------------------------------
   ## RETURN:
	int: -1 if the strig is not an atom element, the index in the periodic table if so.
   -----------------------------------------------------------------------------
*/
int is_valid_element(const char *str, int ignore_case)
{
	if(str == NULL) return -1 ;
	if(strlen(str) <= 0) return -1 ;

	/* Use temporary variable to work on the string */
	int i ;
	char str_tmp[strlen(str)+1] ;
	strcpy(str_tmp, str) ;

	/* Remove spaces and case if asked*/
	str_trim(str_tmp) ;
	if(ignore_case == 1) {
		str_tmp[0] = tolower(str_tmp[0]) ;
		str_tmp[1] = tolower(str_tmp[1]) ;
	}

	/* Loop over */
	for (i = 0; i < ST_nelem ; i++) {
		char tmp[3] ;
		tmp[0] = ST_pte_symbol[i][0] ;
		tmp[1] = ST_pte_symbol[i][1] ;

		/* Remove case if asked */
		if(ignore_case == 1) {
			tmp[0] = tolower(tmp[0]) ;
			tmp[1] = tolower(tmp[1]) ;
		}
		tmp[2] = '\0' ;

		/* Do the comparison*/
		if(strcmp(str_tmp, tmp) == 0) return i ;
	}
	
	return -1 ;
}
