
#include "../headers/rpdb.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					rpdb.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			01-04-08
##
## ----- SPECIFICATIONS
##
##  PDB utilities, mainly reading routines. Not all information are stored,
##  this reader is specific to fpocket, and essentially deals with ATOMS,
##  HETATMS and several other fields.
##
##
##  Curently:
##		- Hydrogens present in the PDB are kept
##		- HETATM are ignored except for specific cofactor, small molecule... 
##		  listed in ST_keep_hetatm variable, and for a  given ligand, defined by 
##		  its resname.
##		- Solvent molecules are ignored
##
##	The reading process stop at the end of the file or as soon as the END flag
##	is reached.
##
##  We used VMD as a base for the code developpement, but now codes are quite
##  differents... 
##
## ----- MODIFICATIONS HISTORY
##
##  17-03-09    (v)  Improved atom type guessing
##  10-03-09    (v)  Atom type guessed using resname when element symbol is missing
##  11-02-09    (v)  Added list of pointer on all atoms (usefull for sorting)
##  22-01-09    (p)  Eliminate double entries in the cofactor list
##	19-01-09	(v)  Open pdb file checking the case
##  28-11-08    (v)  Comments UTD + minor corrections
##  20-11-08    (p)  Adding support of reading PDB multiple occupancies (only 
##					 the first is read) & MSE added to kept HETATM list
##	03-11-08	(v)  Code slightly restructured, comments added
##	01-04-08	(v)  Added template for comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## ----- TODO or SUGGESTIONS
##
##	(v) Improve element guessing by checking resname.
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
	A list of HETATM to keep in the PDB structure.
	REF Peter?
*/

static const char *ST_keep_hetatm[] = {

	"HEA", "HBI", "BIO", "CFM", "CLP", "FES", "F3S", "FS3", "FS4", "BPH",
	"BPB", "BCL", "BCB", "COB",  "ZN", "FEA", "FEO", "H4B", "BH4", "BHS",
	"HBL", "THB", "DDH", "DHE", "HAS", "HDD", "HDM", "HEB", "HEC", "HEO",
	"HES", "HEV", "MHM", "SRM", "VER", "1FH", "2FH", "HC0", "HC1", "HF3",
	"HF5", "NFS", "OMO", "PHF", "SF3", "SF4", "CFM", "CFN", "CLF", "CLP",
	"CN1", "CNB", "CNF", "CUB", "CUM", "CUN", "CUO",  "FS2","FSO", "FSX",
    "PHO", "BH1", "CHL", "CL1", "CL2", "CLA", "CCH", "CFO", "FE2", "FCI", 
    "FCO", "FDC", "FEA", "FEO", "FNE", "HIF", "OFO", "PFC", "HE5", "BAZ", 
    "BOZ",  "FE", "HEM", "HCO", "1CP", "CLN", "COH", "CP3", "DEU", "FDD",
    "FDE", "FEC", "FMI", "HEG", "HNI", "MMP", "MNH", "MNR", "MP1", "PC3",
    "PCU", "PNI", "POR", "PP9", "MSE"
} ;

static const int ST_nb_keep_hetatm = 105 ;

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	rpdb_extract_pdb_atom
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Extract all information given in a pdb ATOM or HETATM line, and store them 
	in given pointers. User must therefore provide enough memory in parameter.
	PDB last known standart:

	COLUMNS      DATA TYPE        FIELD      DEFINITION
	------------------------------------------------------
	1 -  6      Record name      "ATOM    "
	7 - 11      Integer          serial     Atom serial number.
	13 - 16      Atom             name       Atom name.
	17           Character        altLoc     Alternate location indicator.
	18 - 20      Residue name     resName    Residue name.
	22           Character        chainID    Chain identifier.
	23 - 26      Integer          resSeq     Residue sequence number.
	27           AChar            iCode      Code for insertion of residues.
	31 - 38      Real(8.3)        x          Orthogonal coordinates for X in 
											 Angstroms
	39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in 
											 Angstroms
	47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in 
											 Angstroms
	55 - 60      Real(6.2)        occupancy  Occupancy.
	61 - 66      Real(6.2)        tempFactor Temperature factor.
	77 - 78      LString(2)       element    Element symbol, right-justified.
	79 - 80      LString(2)       charge     Charge on the atom.

   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ char *pdb_line	: The PDB line containings info
	@ int *atm_id		: Pointer to atom ID
	@ char *name		: Pointer to atom name
	@ char *res_name	: Pointer to residue name
	@ char *chain		: Pointer to chain name
	@ char *seg_name	: Pointer to segment	
	@ int *res_id 		: Pointer to residue ID
	@ char *insert		: Pointer to insertion code
	@ char *alt_loc		: Pointer to alternate location
	@ char *elem_symbol	: Pointer to element symbol
	@ float *x, *y, *z	: Pointer to coordinates
	@ float *occ		: Pointer to occupency
	@ float *bfactor	: Pointer to b-factor
	@ char *symbol		: Pointer to symbol
	@ float *bfactor	: Pointer to charge
 	@ int guess_flag	: Flag if elements were guessed
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/

void rpdb_extract_pdb_atom( char *pdb_line, char *type, int *atm_id, char *name, 
							char *alt_loc, char *res_name, char *chain, 
							int *res_id, char *insert, 
							float *x, float *y, float *z, float *occ, 
							float *bfactor, char *symbol, int *charge, int *guess_flag)
{
/* Position:          1         2         3         4         5         6 */
/* Position: 123456789012345678901234567890123456789012345678901234567890 */
/* Record:   ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 */

/* Position: 6         7         8 */
/* Position: 012345678901234567890 */
/* Record:   0 11.92           N   */

	int rlen = strlen(pdb_line) ;

	char *prt,
		 ctmp ;

	/* Record type */
	strncpy(type, pdb_line, 6) ;

	/* Atom ID */
	prt = pdb_line + 6 ; 
	ctmp = pdb_line[11] ; pdb_line[11] = '\0' ;
	*atm_id = atoi(prt) ; pdb_line[11] = ctmp ;
	
	/* Atom name */
	strncpy(name, pdb_line + 12, 4);
	name[4] = '\0';
	str_trim(name) ;

	/* Alternate location identifier */
	*alt_loc = pdb_line[16] ;
	
	/* Residue name */
	rpdb_extract_atm_resname(pdb_line, res_name) ;

	/* Chain name */
	chain[0] = pdb_line[21];
	chain[1] = '\0';

	/* Residue id number */
	prt = pdb_line + 22 ; 
	ctmp = pdb_line[26] ; pdb_line[26] = '\0' ;
	*res_id = atoi(prt) ; pdb_line[26] = ctmp ;

	/* Insertion code */
	*insert = pdb_line[26];

	/* x, y, and z coordinates, occupancy and b-factor */
	rpdb_extract_atom_values(pdb_line, x, y, z, occ, bfactor);

	/* Atomic element symbol (if does not exists, guess it based on
	 * atom name */
	if (rlen >= 77) {
		strncpy(symbol, pdb_line + 76, 2);
		symbol[2] = '\0';
		str_trim(symbol); /* remove spaces */
		if(strlen(symbol) < 1) {
			guess_element(name, symbol) ;
                        *guess_flag+=1;
		}
	}
	else {
		guess_element(name, symbol) ;
                *guess_flag+=1;
	}
	str_trim(symbol); /* remove spaces */
	
	/* Charge */
	if(rlen >= 79) {
            char buf[4] = "   " ;
		if((pdb_line[78] == ' ' && pdb_line[79] == ' ') || pdb_line[78] == '\n'){
			*charge = 0 ;
		}
		else {
			buf[0] = pdb_line[78] ;
			buf[1] = pdb_line[79] ;
			buf[2] = '\0' ;
			*charge = (int) atoi(buf) ;
		}
	}
	else *charge = 0 ;
	
}


/**-----------------------------------------------------------------------------
   ## FUNCTION:
	guess_element
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Guess the element of the atom based on atom name. The pattern matched here
	have been taken from the MOE PDB reader.

 	' CL#' => Chlorine (use El field if it matches ' @@?' && HETATAM)
	CH2T from CT !!! => what not

	const AtomPatterns = [
	  N:'[A-G,I-L,N-Z]N#*',
	  O:['[A-B,D-G,I-L,N-Z]O*','OP[A-C]#','CO[A-Z,0-9]*','OE##'],
	  P:'[A-G,I-K,M-N,P-Z]P*',
	  C:['[A-G,I-Z]C#*','C[B-G,I-K,M,P-T,V-Z]#*','#CH#', 'BC  '],
	  H:['H[0-9,A-E,H-Z]*','#[0-9,A-Z]H*','?H[A-Z,0-9]*', 'HG##'],
	  CL:'#CL#',        // check this
	  S:'[P,N]S#*',
	  SE:'NSE1'
	];
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ char *atom_name	: The atom name
	@ char *res_name	: OUTPUT the element guessed
   -----------------------------------------------------------------------------
   ## RETURN:
	void (element is the output)
   -----------------------------------------------------------------------------
*/
void guess_element(char *aname, char *element)
{
	/* Use a temporary variable for atomname, mainly to remove spaces */
	char tmp[strlen(aname)+1] ;
	strcpy(tmp, aname) ;
	
	str_trim(tmp) ;

	char *ptmp = tmp ;

	/* Move to the second caracter if we find a number */ 
        if(isdigit(tmp[0])) {
            ptmp = ptmp+1 ;
        }
        else {
            element[0] = ' ';
            element[1] = ptmp[0] ;
            element[2] = '\0';
            return;
        }
	/* Check if its a valid element */
	int index = is_valid_element(ptmp, 1) ;
	if(index != -1) {
		strcpy(element, ptmp );
		return ;
	}

	/* Here we have a special case... So take the first and second */
	element[0] = ' ';
	element[1] = ptmp[0] ;
	element[2] = '\0';
}

int is_N(char *aname)
{	/* N:'[A-G,I-L,N-Z]N#*' */

	if(aname[0] == 'N' && isdigit(aname[1])) return 1 ;
	if( aname[0] != 'H' && aname[0] != 'M' && aname[1] == 'N'
			&& str_is_number(aname, 0)) return 1 ;

	return 0 ;
}

int is_O(char *aname)
{
/*
	  O:['[A-B,D-G,I-L,N-Z]O*','OP[A-C]#','CO[A-Z,0-9]*','OE##']
*/
	if(aname[0] == 'O') {
		/* TESTING     'OP[A-C]#' */
		if( aname[1] == 'P' && (aname[2] == 'A' || aname[2] == 'B' || aname[2] == 'C')
			&& isdigit(aname[3])) {
			return 1 ;
		}

		/* TESTING     'OE##' */
		if(aname[1] == 'E' && isdigit(aname[2]) && isdigit(aname[3])) return 1 ;
	}
	else {
		/* TESTING     '[A-B,D-G,I-L,N-Z]O*' */
		if( aname[0] != 'C' && aname[0] != 'H' && aname[0] != 'M' && aname[1] == 'O')
			return 1 ;

		/* TESTING     'CO[A-Z,0-9]*' */
		if(aname[0] == 'C' && aname[1] == 'O' && aname[2] != ' ' && aname[3] != ' ')
			return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	rpdb_extract_atm_resname
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Extract the residu name for an ATOM or HETATM pdb record. To remember:

	COLUMNS      DATA TYPE        FIELD      DEFINITION
	------------------------------------------------------
	18 - 20      Residue name     resName    Residue name.

	The memory to store the name has to be provided by the user.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ char *pdb_line	: The PDB line containings info
	@ char *res_name	: Pointer to residue name
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void rpdb_extract_atm_resname(char *pdb_line, char *res_name)
{
/* Position:          1         2         3         4         5         6 */
/* Position: 123456789012345678901234567890123456789012345678901234567890 */
/* Record:   ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 */

/* Position: 6         7         8 */
/* Position: 012345678901234567890 */
/* Record:   0 11.92           N   */
	
/* Residue name */
	strncpy(res_name, pdb_line + 17, 4);
	res_name[4] = '\0';
	/*str_trim(res_name); */
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	rpdb_extract_atom_values
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Extract coordinates, occupancy and bfactor values from a pdb ATOM or HETATM
	line, and store them in given pointers.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ char *pdb_line	: The PDB line containings info
	@ float *x, *y, *z	: Pointer to coordinates
	@ float *occ		: Pointer to occupency
	@ float *bfactor	: Pointer to b-factor
   -----------------------------------------------------------------------------
   ## RETURN: void
   -----------------------------------------------------------------------------
*/
void rpdb_extract_atom_values(char *pdb_line, float *x, float *y, float *z,
							  float *occ, float *bfactor) 
{
/* Position:          1         2         3         4         5         6 */
/* Position: 123456789012345678901234567890123456789012345678901234567890 */
/* Record:   ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 */

/* Position: 6         7         8 */
/* Position: 012345678901234567890 */
/* Record:   0 11.92           N   */
	
	char *ptr,
		 ctmp ;
	
	ptr = pdb_line + 30 ;
	ctmp = pdb_line[38] ; pdb_line[38] = '\0' ;
	*x = (float) atof(ptr) ; pdb_line[38] = ctmp ;

	ptr = pdb_line + 38 ;
	ctmp = pdb_line[46] ; pdb_line[46] = '\0' ;
	*y = (float) atof(ptr) ; pdb_line[46] = ctmp ;

	ptr = pdb_line + 46 ;
	ctmp = pdb_line[54] ; pdb_line[54] = '\0' ;
	*z = (float) atof(ptr) ; pdb_line[54] = ctmp ;

	ptr = pdb_line + 54 ;
	ctmp = pdb_line[60] ; pdb_line[60] = '\0' ;
	*occ = (float) atof(ptr) ; pdb_line[60] = ctmp ;
	
	ptr = pdb_line + 60 ;
	ctmp = pdb_line[66] ; pdb_line[66] = '\0' ;
	*bfactor = (float) atof(ptr) ; pdb_line[66] = ctmp ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	rpdb_extract_cryst1
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Extract information on a box size from a pdb CRYSTL line, and store them 
	in given pointers.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ char *pdb_line				: The PDB line containings info
	@ float *alpha, *beta, *gamma	: Pointer to angles
	@ float *A, B, C				: Pointer sides length
   -----------------------------------------------------------------------------
   ## RETURN: void
   -----------------------------------------------------------------------------
*/
void rpdb_extract_cryst1(char *pdb_line, float *alpha, float *beta, float *gamma, 
						 float *a, float *b, float *c) 
{
/* Position:          1         2         3         4         5         6 */
/* Position: 123456789012345678901234567890123456789012345678901234567890 */
/* Record:   ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 */

/* Position: 6         7         8 */
/* Position: 012345678901234567890 */
/* Record:   0 11.92           N   */

	char ch, *s;
	
	s = pdb_line+6 ;
	ch = pdb_line[15] ; pdb_line[15] = '\0' ;
	*a = (float) atof(s) ;

	s = pdb_line+15 ;
	*s = ch ; ch = pdb_line[24]; pdb_line[24] = '\0' ;
	*b = (float) atof(s) ;

	s = pdb_line+24 ;
	*s = ch; ch = pdb_line[33]; pdb_line[33] = '\0' ;
	*c = (float) atof(s) ;

	s = pdb_line+33;
	*s = ch; ch = pdb_line[40]; pdb_line[40] = '\0' ;
	*alpha = (float) atof(s) ;

	s = pdb_line+40;
	*s = ch; ch = pdb_line[47]; pdb_line[47] = '\0' ;
	*beta = (float) atof(s) ;

	s = pdb_line+47;
	*s = ch; ch = pdb_line[54]; pdb_line[54] = '\0' ;
	*gamma = (float) atof(s) ;
}


/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	rpdb_open
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Open a PDB file, alloc memory for all information on this pdb, and store 
	several information like the number of atoms, the header, the remark... 
	This first reading of PDB rewinds the FILE* pointer. No coordinates are
	actually read.

	Hydrogens are conserved.
	All HETATM are removed, except the given ligand if we have to keep it, and
	important HETATM listed in the static structure at the top of this file.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ const char *fpath    : The pdb path.
	@ const char *ligan    : Ligand resname.
	@ const char *keep_lig :  Keep the given ligand or not?
   -----------------------------------------------------------------------------
   ## RETURN: 
	s_pdb: data containing PDB info.
   -----------------------------------------------------------------------------
*/
s_pdb* rpdb_open(char *fpath, const char *ligan, const int keep_lig)
{
	s_pdb *pdb = NULL ;

	char buf[M_PDB_BUF_LEN],
		 resb[5] ;
	
	int nhetatm = 0,
		natoms = 0,
		natm_lig = 0 ;
	int i ;
	
	pdb = (s_pdb *) my_malloc(sizeof(s_pdb)) ; ;
	
	/* Open the PDB file in read-only mode */
	pdb->fpdb = fopen_pdb_check_case(fpath, "r");
	if (!pdb->fpdb) {
		my_free(pdb) ;
		fprintf(stderr, "! File %s does not exist\n", fpath) ;
		return NULL ;
	}

	while(fgets(buf, M_PDB_LINE_LEN + 2, pdb->fpdb)) {
		if (!strncmp(buf, "ATOM ",  5)) {
		/* Check if this is the first occurence of this atom*/
			if(buf[16]==' ' || buf[16]=='A'){  
			/* Atom entry: check if there is a ligand in there (just in case)... */
				rpdb_extract_atm_resname(buf, resb) ;
				if( ligan && ligan[0] == resb[0] && ligan[1] == resb[1] 
					&& ligan[2] == resb[2]){

					if(keep_lig) {
						natm_lig ++ ;
						natoms++ ;
					}
				}
				else {
					natoms++ ;
				}
			}
		}
		else if(!strncmp(buf, "HETATM", 6)) {
			/*Check again for the first occurence*/
			if(buf[16]==' ' || buf[16]=='A'){ 
				/* Hetatom entry: check if there is a ligand in there too... */
				rpdb_extract_atm_resname(buf, resb) ;
				if( keep_lig && ligan && ligan[0] == resb[0] && ligan[1] == resb[1] 
					&& ligan[2] == resb[2]){
					natm_lig ++ ; natoms++ ;
				}
				else {
				/* Keep specific HETATM given in the static list ST_keep_hetatm */
					for(i = 0 ; i < ST_nb_keep_hetatm ; i++) {
							if(ST_keep_hetatm[i][0] == resb[0] && ST_keep_hetatm[i][1]
								== resb[1] && ST_keep_hetatm[i][2] == resb[2]) {
									nhetatm++ ; natoms++ ;
									break ;
							}
					}
				}
			}
		}
/*
		else if (!strncmp(buf, "HEADER", 6)) 
			strncpy(pdb->header, buf, M_PDB_BUF_LEN) ;
*/
		
		else if (!strncmp(buf, "END", 3)) break ;
	}

	if (natoms == 0) {
		fprintf(stderr, "! File '%s' contains no atoms...\n", fpath) ;
		my_free(pdb) ;
	
		return NULL ;
	}

	/* Alloc needed memory */
	pdb->latoms = (s_atm*) my_calloc(natoms, sizeof(s_atm)) ;
	pdb->latoms_p = (s_atm**) my_calloc(natoms, sizeof(s_atm*)) ;

	if(nhetatm > 0) pdb->lhetatm = (s_atm**) my_calloc(nhetatm, sizeof(s_atm*)) ;
	else pdb->lhetatm = NULL ;
	
	if(natm_lig > 0) pdb->latm_lig = (s_atm**) my_calloc(natm_lig, sizeof(s_atm*)) ;
	else pdb->latm_lig = NULL ;
	
	pdb->natoms = natoms ;
	pdb->nhetatm = nhetatm ;
	pdb->natm_lig = natm_lig ;
	rewind(pdb->fpdb) ;

	return pdb ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	rpdb_read
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Read and store information on atoms for a pdb file.
    Curently:
		- Hydrogens present in the PDB are kept
		- HETATM are ignored except for specific cofactor, small molecule... 
		  listed in ST_keep_hetatm variable, and for a  given ligand, defined by 
		  its resname.
		- Solvent molecules are ignored
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pdb *pdb           : The structure to fill
	@ const char *ligand   : The ligand resname
	@ const char *keep_lig :  Keep the given ligand or not?
   -----------------------------------------------------------------------------
   ## RETURN:
   -----------------------------------------------------------------------------
*/
void rpdb_read(s_pdb *pdb, const char *ligan, const int keep_lig) 
{
	int i,
		iatoms,
		ihetatm,
		iatm_lig,
		ligfound ;

	char pdb_line[M_PDB_BUF_LEN],
		 resb[5] ;					/* Buffer for the current residue name */

	s_atm *atom = NULL ;
	s_atm *atoms = pdb->latoms ;
	s_atm **atoms_p = pdb->latoms_p ;
	s_atm **atm_lig = pdb->latm_lig ;
        int guess_flag=0;
	iatoms = 0 ;
	ihetatm = 0 ;
	iatm_lig = 0 ;
	ligfound = 0 ;

	/* Loop over the pdb file */ 
	while(fgets(pdb_line, M_PDB_LINE_LEN + 2, pdb->fpdb)) {
		if (strncmp(pdb_line, "ATOM ",  5) == 0) {
			if(pdb_line[16]==' ' || pdb_line[16]=='A'){ /*if within first occurence*/
				/* Store ATOM entry */
				rpdb_extract_atm_resname(pdb_line, resb) ;
				/* Check if the desired ligand is in such entry */
				if( ligan && ligan[0] == resb[0] && ligan[1] == resb[1] 
					&& ligan[2] == resb[2]){
					if(keep_lig) {
						atom = atoms + iatoms ;
						
						/* Read atom information */
						rpdb_extract_pdb_atom(pdb_line, atom->type, &(atom->id), 
								atom->name, &(atom->pdb_aloc), atom->res_name, 
								atom->chain, &(atom->res_id), &(atom->pdb_insert), 
								&(atom->x), &(atom->y), &(atom->z), 
								&(atom->occupancy), &(atom->bfactor), atom->symbol,
								&(atom->charge), &guess_flag);

						/* Store additional information not given in the pdb */
						atom->mass = pte_get_mass(atom->symbol) ;
						atom->radius = pte_get_vdw_ray(atom->symbol) ;
						atom->electroneg = pte_get_enegativity(atom->symbol) ;
						atom->sort_x = -1 ;
						
						atoms_p[iatoms] = atom ;
						iatoms++ ;

						atm_lig[iatm_lig] = atom ;
						iatm_lig ++ ;
						ligfound = 1 ;
					}
				}
				else {
				/* A simple atom not supposed to be stored as a ligand */
					atom = atoms + iatoms ;
					rpdb_extract_pdb_atom(pdb_line, atom->type, &(atom->id), 
							atom->name, &(atom->pdb_aloc), atom->res_name, 
							atom->chain, &(atom->res_id), &(atom->pdb_insert), 
							&(atom->x), &(atom->y), &(atom->z), &(atom->occupancy), 
							&(atom->bfactor), atom->symbol, &(atom->charge), &guess_flag);

					/* Store additional information not given in the pdb */
					atom->mass = pte_get_mass(atom->symbol) ;
					atom->radius = pte_get_vdw_ray(atom->symbol) ;
					atom->electroneg = pte_get_enegativity(atom->symbol) ;
					atom->sort_x = -1 ;

					atoms_p[iatoms] = atom ;
					iatoms++ ;
				}
			}
		}
		else if(strncmp(pdb_line, "HETATM", 6) == 0) {
			if(pdb_line[16]==' ' || pdb_line[16]=='A'){/*first occurence*/
				/* Check HETATM entry */
				rpdb_extract_atm_resname(pdb_line, resb) ;
				/* Check if the desired ligand is in HETATM entry */
				if( ligan && keep_lig && ligan[0] == resb[0] && ligan[1] == resb[1] 
					&& ligan[2] == resb[2]){

					atom = atoms + iatoms ;
					rpdb_extract_pdb_atom(pdb_line, atom->type, &(atom->id), 
							atom->name, &(atom->pdb_aloc), atom->res_name, 
							atom->chain, &(atom->res_id), &(atom->pdb_insert), 
							&(atom->x), &(atom->y), &(atom->z), &(atom->occupancy), 
							&(atom->bfactor), atom->symbol, &(atom->charge), &guess_flag);

					/* Store additional information not given in the pdb */
					atom->mass = pte_get_mass(atom->symbol) ;
					atom->radius = pte_get_vdw_ray(atom->symbol) ;
					atom->electroneg = pte_get_enegativity(atom->symbol) ;
					atom->sort_x = -1 ;

					atoms_p[iatoms] = atom ;
					atm_lig[iatm_lig] = atom ;
					
					iatm_lig ++ ; iatoms++ ;
					ligfound = 1 ;
				}
				else if(pdb->lhetatm) {
				/* Keep specific HETATM given in the static list ST_keep_hetatm. */
					for(i = 0 ; i < ST_nb_keep_hetatm ; i++) {
						if( ST_keep_hetatm[i][0] == resb[0] && ST_keep_hetatm[i][1] 
							== resb[1] && ST_keep_hetatm[i][2] == resb[2]) {
							atom = atoms + iatoms ;
							rpdb_extract_pdb_atom(pdb_line, atom->type, &(atom->id), 
									atom->name, &(atom->pdb_aloc), atom->res_name, 
									atom->chain, &(atom->res_id), &(atom->pdb_insert), 
									&(atom->x), &(atom->y), &(atom->z), 
									&(atom->occupancy), &(atom->bfactor), 
									atom->symbol, &(atom->charge), &guess_flag);

						/* Store additional information not given in the pdb */
							atom->mass = pte_get_mass(atom->symbol) ;
							atom->radius = pte_get_vdw_ray(atom->symbol) ;
							atom->electroneg = pte_get_enegativity(atom->symbol) ;
							atom->sort_x = -1 ;
							
							atoms_p[iatoms] = atom ;
							pdb->lhetatm[ihetatm] = atom ;
							ihetatm ++ ; iatoms++ ;
							break ;
						}
					}
				}
			}
		}
		else if (strncmp(pdb_line, "CRYST1",  6) == 0)  {
			rpdb_extract_cryst1(pdb_line, &(pdb->alpha), &(pdb->beta), &(pdb->gamma),
									 &(pdb->A),  &(pdb->B), &(pdb->C));
		}
		else if (!strncmp(pdb_line, "END", 3)) break ;
	}

       
        if(guess_flag>0) {
            fprintf(stderr, ">! Warning: You did not provide a standard PDB file.\nElements were guessed by fpocket, because not provided in the PDB file. \nThere is no guarantee on the results!\n");
        }
        
	if(ligan && keep_lig && (ligfound == 0 || pdb->natm_lig <= 0)) {
		fprintf(stderr, ">! Warning: ligand '%s' not found in the pdb...\n", ligan) ;
		if(pdb->latm_lig) fprintf(stderr, "! Ligand list is not NULL however...\n") ;
		if(ligfound == 1) fprintf(stderr, "! And ligfound == 1!! :-/\n") ;
	}
	else if(ligfound == 1 && iatm_lig <= 0) {
		fprintf(stderr, ">! Warning: ligand '%s' has been detected but no atoms \
						has been stored!\n", ligan) ;
	}
	else if((ligfound == 1 && pdb->natm_lig <= 0) || (pdb->natm_lig <=0 
			 && iatm_lig > 0)) {
		fprintf(stderr, ">! Warning: ligand '%s' has been detected in rpdb_read \
						but not in rpdb_open!\n", ligan) ;
	}

}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	free_pdb_atoms
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Free memory for s_pdb structure
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pdb *pdb: pdb struct to free
   -----------------------------------------------------------------------------
   ## RETURN:
	void
   -----------------------------------------------------------------------------
*/
void free_pdb_atoms(s_pdb *pdb) 
{
	if(pdb) {
		if(pdb->lhetatm) {
			my_free(pdb->lhetatm) ;
			pdb->lhetatm = NULL ;		
		}
		if(pdb->latoms) {
			my_free(pdb->latoms) ;
			pdb->latoms = NULL ;
		}
		if(pdb->latm_lig) {
			my_free(pdb->latm_lig) ; 
			pdb->latm_lig = NULL ;
		}
		if(pdb->fpdb) {
			fclose(pdb->fpdb) ;
			pdb->fpdb = NULL ;	
		}
		
		if(pdb->latoms_p) {
			my_free(pdb->latoms_p) ;
			pdb->latoms_p = NULL ;
		}

		my_free(pdb) ;
	}
}
