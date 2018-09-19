 
#include "../headers/fparams.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					fparams.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			28-11-08
##
## ----- SPECIFICATIONS
##
##	Handle parameters (parse the command line and sore values)
##	for the fpocket programm.
##
## ----- MODIFICATIONS HISTORY
##
-##	17-03-09	(v)  Segfault avoided when freeing pdb list
##	15-12-08	(v)  Added function to check if a single letter is a fpocket
##					 command line option (usefull for t/dpocket) + minor modifs
##	28-11-08	(v)  List of pdb taken into account as a single file input.
##					 Comments UTD
##	27-11-08	(v)  PDB file check moved in fpmain + minor modif + relooking
##	01-04-08	(v)  Added comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## ----- TODO or SUGGESTIONS
##
##	(v) Check and update if necessary comments of each function!!
##	(v) Review the main function and handle all possible crashes.
##

*/

/**-----------------------------------------------------------------------------
   ## FUNCTION:
	init_def_fparams
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Initialisation of default parameters
   -----------------------------------------------------------------------------
   ## PARAMETRES: void
   -----------------------------------------------------------------------------
   ## RETURN: 
	s_fparams*: Pointer to allocated paramers.
   -----------------------------------------------------------------------------
*/
s_fparams* init_def_fparams(void)
{
	s_fparams *par = (s_fparams *) my_malloc(sizeof(s_fparams)) ;

	par->min_apol_neigh = M__MIN_APOL_NEIGH_DEFAULT ;
	par->asph_min_size = M_MIN_ASHAPE_SIZE_DEFAULT ;
	par->asph_max_size = M_MAX_ASHAPE_SIZE_DEFAULT ;
	par->sl_clust_max_dist = M_SLCLUST_MAX_DIST ;
	par->sl_clust_min_nneigh = M_SLCLUST_MIN_NUM_NEIGH ;
	par->pdb_path[0] = 0 ;
	par->basic_volume_div = M_BASIC_VOL_DIVISION ;
	par->nb_mcv_iter = M_MC_ITER ;
	par->min_pock_nb_asph = M_MIN_POCK_NB_ASPH ;
	par->refine_clust_dist = M_REFINE_DIST ;
	par->refine_min_apolar_asphere_prop = M_REFINE_MIN_PROP_APOL_AS ;
	par->clust_max_dist = M_CLUST_MAX_DIST ;
	par->npdb = 0 ;
	par->pdb_lst = NULL ;

	return par ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	get_fpocket_args
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	This function analyse the user's command line and parse it to store parameters
	for the pocket finder programm.
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ int nargs   :  Number of arguments
	@ char **args : Arguments of main program
   -----------------------------------------------------------------------------
   ## RETURN: 
	s_params*: Pointer to parameters
   -----------------------------------------------------------------------------
*/
s_fparams* get_fpocket_args(int nargs, char **args)
{
	int i,
		npdb = 0,
		status = 0 ;

	s_fparams *par = init_def_fparams() ;
	char *pdb_lst = NULL ;
	
	//read arguments by flags
	for (i = 1; i < nargs; i++) {
		if (strlen(args[i]) == 2 && args[i][0] == '-' && i < (nargs-1)) {
			switch (args[i][1]) {
				case M_PAR_MAX_ASHAPE_SIZE	  : 
					status += parse_asph_max_size(args[++i], par) ;		break ;
				case M_PAR_MIN_ASHAPE_SIZE	  : 
					status += parse_asph_min_size(args[++i], par) ;		break ;
				case M_PAR_MIN_APOL_NEIGH	  : 
					status += parse_min_apol_neigh(args[++i], par) ;	break ;
				case M_PAR_CLUST_MAX_DIST	  : 
					status += parse_clust_max_dist(args[++i], par) ;	break ;
				case M_PAR_SL_MAX_DIST		  : 
					status += parse_sclust_max_dist(args[++i], par) ;	break ;
				case M_PAR_SL_MIN_NUM_NEIGH   : 
					status += parse_sclust_min_nneigh(args[++i], par) ; break ;
				case M_PAR_MC_ITER 			  : 
					status += parse_mc_niter(args[++i], par) ;			break ;
				case M_PAR_BASIC_VOL_DIVISION : 
					status += parse_basic_vol_div(args[++i], par) ;		break ;
				case M_PAR_MIN_POCK_NB_ASPH   : 
					status += parse_min_pock_nb_asph(args[++i], par) ;	break ;
				case M_PAR_REFINE_DIST		  : 
					status += parse_refine_dist(args[++i], par) ;		break ;
				case M_PAR_REFINE_MIN_NAPOL_AS: 
					status += parse_refine_minaap(args[++i], par) ; 
					break ;
				case M_PAR_PDB_LIST :
					pdb_lst = args[++i] ; break ;
					
				case M_PAR_PDB_FILE			  : 
						if(npdb >= 1) fprintf(stderr, 
							"! Only first input pdb will be used.\n") ;
						else {
							strcpy(par->pdb_path, args[++i]) ; npdb++ ;
						}
						break ;
				default: break ;
			}
		}
	}
	
	if(status > 0) {
 		free_fparams(par) ;
		print_pocket_usage(stdout);
		return NULL ;
	}
	
	par->npdb = npdb ;
	
	/* Handle a file containing a list of PDB */
	if(pdb_lst != NULL) {
		FILE *f = fopen(pdb_lst, "r") ;
		if(f != NULL) {
			/* Count the number of lines */
			int n = 0 ;
			char cline [M_MAX_PDB_NAME_LEN + 1] ;
			
			while(fgets(cline, M_MAX_PDB_NAME_LEN, f) != NULL) {
				if(strcmp("\n", cline) != 0) {
					n ++ ;
				}
			}
			fclose(f) ;
			if(n == 0) {
				return par ;
			}
			
			/* Allocate memory and store each line */
			par->pdb_lst = (char **)my_malloc(n*sizeof(char*)) ;

			f = fopen (pdb_lst, "r") ; 
			int i = 0, l = 0 ;
			while(fgets(cline, M_MAX_PDB_NAME_LEN, f) != NULL) {
				if(strcmp("\n", cline) != 0) {
					l = strlen(cline) ;
					if(cline[l-1] == '\n') {
						l -- ;
						cline[l] = '\0' ;
					}
					char *line = (char *) my_malloc((l+1)*sizeof(char)) ;
					memcpy (line, cline, l+1);
	
					par->pdb_lst[i] = line ;
					i ++ ;
				}
			}
			
			par->npdb = n ;
		}
	}
	
	return par;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_clust_max_dist
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the distance criteria first clustering algorithm.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid float), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_clust_max_dist(char *str, s_fparams *p) 
{
	if(str_is_float(str, M_NO_SIGN)) {
		p->clust_max_dist = atof(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the single linkage max dist.\n", str) ;
		return 1 ;
	}

	return 0 ;
}


/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_sclust_max_dist
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the distance criteria in the single linkage clustering.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid float), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_sclust_max_dist(char *str, s_fparams *p) 
{
	if(str_is_float(str, M_NO_SIGN)) {
		p->sl_clust_max_dist = atof(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the single linkage max dist.\n", str) ;
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_sclust_min_nneigh
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the number of neighbours in the single linkage clustering.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid int), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_sclust_min_nneigh(char *str, s_fparams *p) 
{
	if(str_is_number(str, M_NO_SIGN)) {
		p->sl_clust_min_nneigh = atoi(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the single linkage max dist.\n", str) ;
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_min_apol_neigh
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the minimum number of apolar contacted atom for an alpha
	sphere to be considered as apolar.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid int), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_min_apol_neigh(char *str, s_fparams *p)
{
	if(str_is_number(str, M_NO_SIGN)) {
		p->min_apol_neigh = (int) atoi(str) ;
		if(p->min_apol_neigh < 0) p->min_apol_neigh = 0 ;
		if(p->min_apol_neigh > 4) p->min_apol_neigh = 4 ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the min radius of alpha shperes.\n", str) ;
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_asph_min_size
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the minimum radius of each alpha shpere
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str: The string to parse
	@ s_fparams *p: The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid float), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_asph_min_size(char *str, s_fparams *p) 
{
	if(str_is_float(str, M_NO_SIGN)) {
		p->asph_min_size = (float) atof(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the min radius of alpha shperes.\n", str) ;
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_asph_max_size
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the maximum radius of each alpha shpere
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid float), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_asph_max_size(char *str, s_fparams *p) 
{
	if(str_is_float(str, M_NO_SIGN)) {
		p->asph_max_size = (float) atof(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the max radius of alpha shperes.\n", str) ;
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_mc_niter
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the number of iteration for the Monte Carlo volume
	calculation.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid float), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_mc_niter(char *str, s_fparams *p)
{
	if(str_is_float(str, M_NO_SIGN)) {
		p->nb_mcv_iter = (int) atoi(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the number of monte-carlo iteration for the volume.\n", str) ;
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_basic_vol_div
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the number of iteration for the basic volume calculation.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid integer), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_basic_vol_div(char *str, s_fparams *p) 
{
	if(str_is_number(str, M_NO_SIGN)) {
		p->basic_volume_div = (int) atoi(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the precision of the basic volume calculation.\n", str) ;
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_refine_dist
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the distance in the refine algorithm
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid float), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_refine_dist(char *str, s_fparams *p) 
{
	if(str_is_float(str, M_NO_SIGN)) {
		p->refine_clust_dist = (float) atof(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the refine distance.\n", str) ;
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_refine_min_apol
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the minimum number of apolar sphere per pocket.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid integer), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_refine_minaap(char *str, s_fparams *p) 
{
	if(str_is_float(str, M_NO_SIGN)) {
		p->refine_min_apolar_asphere_prop = (float) atof(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the refine distance.\n", str) ;
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	parse_min_pock_nb_asph
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 	
	Parsing function for the minimum number of alpha sphere per pocket.
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *str    : The string to parse
	@ s_fparams *p : The structure than will contain the parsed parameter
   -----------------------------------------------------------------------------
   ## RETURN: 
	int: 0 if the parameter is valid (here a valid integer), 1 if not
   -----------------------------------------------------------------------------
*/
int parse_min_pock_nb_asph(char *str, s_fparams *p) 
{
	if(str_is_number(str, M_NO_SIGN)) {
		p->min_pock_nb_asph = (int) atoi(str) ;
	}
	else {
		fprintf(stdout, "! Invalid value (%s) given for the refine distance.\n", str) ;
		return 1 ;
	}

	return 0 ;
}


/**-----------------------------------------------------------------------------
   ## FUNCTION:
	is_fpocket_opt
   -----------------------------------------------------------------------------
   ## SPECIFICATION:
	Say either or not a single letter code is a fpocket option (excluding
	input file/list option.)
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ const char opt: The one letter code option.
   -----------------------------------------------------------------------------
   ## RETURN:
	integer: 1 if it's a valid option parmeter, 0 if not.
   -----------------------------------------------------------------------------
*/

int is_fpocket_opt(const char opt)
{
	if( opt == M_PAR_MAX_ASHAPE_SIZE ||
		opt == M_PAR_MIN_ASHAPE_SIZE ||
		opt == M_PAR_MIN_APOL_NEIGH ||
		opt == M_PAR_CLUST_MAX_DIST ||
		opt == M_PAR_SL_MAX_DIST ||
		opt == M_PAR_SL_MIN_NUM_NEIGH ||
		opt == M_PAR_MC_ITER ||
		opt == M_PAR_BASIC_VOL_DIVISION ||
		opt == M_PAR_MIN_POCK_NB_ASPH ||
		opt == M_PAR_REFINE_DIST ||
		opt == M_PAR_REFINE_MIN_NAPOL_AS) {
		return 1 ;
	}

	return 0 ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	free_fparams
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Free parameters
   -----------------------------------------------------------------------------
   ## PARAMETRES: 
	@ s_params *p: Pointer to the structure to free
   -----------------------------------------------------------------------------
   ## RETURN: 
	void
   -----------------------------------------------------------------------------
*/
void free_fparams(s_fparams *p)
{
	if(p) {
		if(p->npdb > 0 && p->pdb_lst != NULL) {
			int i ;
			for (i = 0 ; i < p->npdb ; i++) {
				if(p->pdb_lst[i] != NULL) my_free(p->pdb_lst[i]) ;
			}
			my_free(p->pdb_lst) ;
		}
		my_free(p) ;
	}
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	print_pocket_usage
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Displaying usage of the programm in the given buffer
   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ FILE *f: buffer to print in
   -----------------------------------------------------------------------------
   ## RETURN:
    void
   -----------------------------------------------------------------------------
*/
void print_pocket_usage(FILE *f)
{
	f = (f == NULL) ? stdout:f ;

	fprintf(f, M_FP_USAGE) ;
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	print_fparams
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Print function
   -----------------------------------------------------------------------------
   ## PARAMETRES:
    @ s_fparams *p : Parameters to print
	@ FILE *f      : Buffer to write in
   -----------------------------------------------------------------------------
   ## RETURN: 
   -----------------------------------------------------------------------------
*/
void print_fparams(s_fparams *p, FILE *f)
{
	if(p) {
		fprintf(f, "==============\nParameters of fpocket: \n");
		fprintf(f, "> Minimum alpha sphere radius: %f\n", p->asph_min_size);
		fprintf(f, "> Maximum alpha sphere radius: %f\n", p->asph_max_size);
		fprintf(f, "> Minimum number of apolar neighbor: %d\n", p->min_apol_neigh);
		fprintf(f, "> Maximum distance for first clustering algorithm: %f \n", p->clust_max_dist) ;
		fprintf(f, "> Single linkage clustering distance: %f\n", p->sl_clust_max_dist);
		fprintf(f, "> Single linkage clustering neighbor: %d\n", p->sl_clust_min_nneigh);
		fprintf(f, "> Refine clustering distance: %f\n", p->refine_clust_dist);
		fprintf(f, "> Min number of apolar sphere in refine to keep a pocket: %f\n", p->refine_min_apolar_asphere_prop) ;
		fprintf(f, "> Monte carlo iterations: %d\n", p->nb_mcv_iter);
		fprintf(f, "> Basic method for volume calculation: %d\n", p->basic_volume_div);
		fprintf(f, "> PDB file: %s\n", p->pdb_path);
		fprintf(f, "==============\n");
	}
	else fprintf(f, "> No parameters detected\n");
}

