
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

#ifndef DH_MEMHANDLER
#define DH_MEMHANDLER

/* -----------------------------INCLUDES--------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>


/* -------------------------------MACROS--------------------------------------*/

#define M_EXIT 1
#define M_CONTINUE 0

/* ------------------------- PUBLIC STRUCTURES ------------------------------ */

/* ----------------------------PROTOTYPES-------------------------------------*/

void* my_malloc(size_t nb) ;
void* my_calloc(size_t nb, size_t s) ;
void* my_realloc(void *ptr, size_t nb) ;
void my_free(void *bloc) ;
void my_exit(void) ;

void free_all(void) ;
void print_ptr_lst(void) ;

#endif
