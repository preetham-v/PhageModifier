
#include "../headers/pscoring.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					pscoring.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			28-11-08
##
## ----- SPECIFICATIONS
##
## This file stores scoring functions for pockets.
##
## ----- MODIFICATIONS HISTORY
##
##	21-01-09	(v) Added new scoring function
##	28-11-08	(v) Created + Comments UTD
##	
## ----- TODO or SUGGESTIONS
##
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
	score_pocket2
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Set a score to a given pocket. The current scoring function has been determined
	using a logistic regression based on an analysis of pocket descriptors.

   -----------------------------------------------------------------------------
   ## PARAMETRES:
	@ s_pocket *pocket: The pocket
   -----------------------------------------------------------------------------
   ## RETURN:
	float: The score
   -----------------------------------------------------------------------------
*/
float score_pocket(s_desc *pdesc)
{
	float score ;

	/**
	 * Data to use for mean-center normalization step: N = 2
	 *                     MEAN     SD
	 * nas_norm            0.183    0.243
	 * apol_asprop_norm    0.402    0.259
	 * mean_loc_hd_norm    0.334    0.261
	 * polarity_score      7.193    4.197
	 * polarity_score_norm 0.291    0.242
	 * as_density          5.194    1.707
	 * as_density_norm     0.387    0.252
	 *
	 */

/*
	Using m 3.0 M 6.0 D 1.73 i 25 we have for the training set this PLS model
	having 4 components

	CURRENT !!!!!!!!!!!!!!!!! SCORING 1
*/
/*
	Perf:
	Scoring function 1:
			  CPP     OVL
	Data      T1/T3 | T1/T3
	------------------------
	Train   : 62/86 - 65/89
	PP holo : 79/92 - 79/90
	PP apo  : 69/90
	Cheng   : 70/85 - 70/100
	Gold    : 69/90 - 71/90
*/


	score =
        -1.50335
       +30.27950 * (float)pdesc->nas_norm
        -3.40435 * (float)pdesc->prop_asapol_norm
       +11.04704 * (float)pdesc->mean_loc_hyd_dens_norm
        +1.18610 * (float)pdesc->polarity_score
        -2.01214 * (float)pdesc->as_density ;

	score =
        -0.65784
       +29.78270 * (float)pdesc->nas_norm
        -4.06632 * (float)pdesc->prop_asapol_norm
       +11.72346 * (float)pdesc->mean_loc_hyd_dens_norm
        +1.16349 * (float)pdesc->polarity_score
        -2.06835 * (float)pdesc->as_density ;


/*
	Using m 3.0 M 6.0 D 1mean_loc_hyd_dens_norm.73 i 25 n 2 we have for the training set this PLS model
	having 4 components

	 SCORING 2
*/
/*
	Perf:
	Scoring function 1:
			  CPP     OVL
	Data      T1/T3 | T1/T3
	------------------------
	Train   : 62/86 - 65/89
	PP holo : 79/92 - 79/90
	PP apo  : 69/90
	Cheng   : 70/85 - 70/100
	Gold    : 69/91 - 71/90
*/

/*
	score =
        -0.65784
       +29.78270 * (float)pdesc->nas_norm
        -4.06632 * (float)pdesc->prop_asapol_norm
       +11.72346 * (float)pdesc->mean_loc_hyd_dens_norm
        +1.16349 * (float)pdesc->polarity_score
        -2.06835 * (float)pdesc->as_density ;
*/
/*
	Using m 3.0 M 6.0 D 1mean_loc_hyd_dens_norm.73 i 25 n 3 we have for the training set this PLS model
	having 4 components

	 SCORING 3
*/
/*
	Perf:
	Scoring function 1:
			  CPP     OVL
	Data      T1/T3 | T1/T3
	------------------------
	Train   : 59/84 - 64/89
	PP holo : 79/94 - 81/94
	PP apo  : 69/90
	Cheng   : 70/85 - 75/100
	Gold    : 71/91 - 72/89
*/

/*
	score =
        -1.48906
       +29.54059 * (float)pdesc->nas_norm
       +10.73666 * (float)pdesc->mean_loc_hyd_dens_norm
        -3.30562 * (float)pdesc->prop_asapol_norm
        +1.15711 * (float)pdesc->polarity_score
        -1.94912 * (float)pdesc->as_density ;
*/

/*	ON GOLD
	Using m 3.0 M 6.0 D 1mean_loc_hyd_dens_norm.73 i 25 n 2 we have for the training set this PLS model
	having 4 components

	 SCORING 4
*/
/*
	Perf:
	Scoring function 1:
			  CPP     OVL
	Data      T1/T3 | T1/T3
	------------------------
	Train   : 59/84 - 64/89
	PP holo : 79/94 - 81/94
	PP apo  : 71/90
	Cheng   : 70/85 - 75/100
	Gold    : 69/91 - 71/90
*/

/*
	score =
        -1.29456
       +33.45117 * (float)pdesc->nas_norm
       +17.78868 * (float)pdesc->mean_loc_hyd_dens_norm
        -5.23046 * (float)pdesc->prop_asapol_norm
        +1.07977 * (float)pdesc->polarity_score
        -2.00073 * (float)pdesc->as_density ;
*/

/*	ON GOLD
	Using m 3.0 M 6.0 D 1mean_loc_hyd_dens_norm.73 i 25 n 3 we have for the training set this PLS model
	having 4 components

	 SCORING 5
*/

/*
	Perf:
	Scoring function 1:
			  CPP     OVL
	Data      T1/T3 | T1/T3
	------------------------
	Train   : 62/86 - 65/89
	PP holo : 79/90 - 79/88
	PP apo  : 67/90
	Cheng   : 70/85 - 70/100
	Gold    : 70/91 - 71/90
*/
/*

	score =
         -2.29256
       +33.86433 * (float)pdesc->nas_norm
       +17.55332 * (float)pdesc->mean_loc_hyd_dens_norm
        -4.90910 * (float)pdesc->prop_asapol_norm
        +1.11252 * (float)pdesc->polarity_score
        -1.88681 * (float)pdesc->as_density;
*/

/*
	score =
        -0.04719
       +27.28918 * (float)pdesc->nas_norm
        -3.28306 * (float)pdesc->prop_asapol_norm
       +11.24130 * (float)pdesc->mean_loc_hyd_dens_norm
        +1.24804 * (float)pdesc->polarity_score
        -2.63044 * (float)pdesc->as_density
        +5.42051 * (float)pdesc->as_max_dst_norm ;
*/



	return score ;
}
