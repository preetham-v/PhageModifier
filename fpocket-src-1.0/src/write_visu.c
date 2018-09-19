#include "../headers/write_visu.h"

/**

## ----- GENERAL INFORMATION
##
## FILE 					write_visu.c
## AUTHORS					P. Schmidtke and V. Le Guilloux
## LAST MODIFIED			02-12-08
##
## ----- SPECIFICATIONS
##
##		Write output script to launch visualization of fpocket output using
##		pymol and VMD.
##
## ----- MODIFICATIONS HISTORY
##
##  02-12-08    (v)  Comments UTD
##  29-11-08    (p)  Enhanced VMD output, corrected bug in pymol output
##  20-11-08    (p)  Just got rid of a memory issue (fflush after fclose) 
##	01-04-08	(v)  Added template for comments and creation of history
##	01-01-08	(vp) Created (random date...)
##	
## ----- TODO or SUGGESTIONS
##
## (v) Handle system command failure
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
	write_visualization
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Write visualization scripts for VMD and PyMol
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *pdb_name     : pdb code
	@ char *pdb_out_name : name of the pdb output file
   -----------------------------------------------------------------------------
   ## RETURN: void
   -----------------------------------------------------------------------------
*/

void write_visualization(char *pdb_name,char *pdb_out_name)
{
	write_vmd(pdb_name,pdb_out_name);
 	write_pymol(pdb_name,pdb_out_name);
}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	write_vmd
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	Write visualization script for VMD
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *pdb_name : pdb code
	@ char *pdb_out_name : name of the pdb output file
   -----------------------------------------------------------------------------
   ## RETURN: void
   -----------------------------------------------------------------------------
*/


void write_vmd(char *pdb_name,char *pdb_out_name)
{
	char fout[250] = "" ;
	char fout2[250] = "" ;
	char sys_cmd[250] ="";
	char c_tmp[255];
	int status ;
	
	strcpy(c_tmp,pdb_name);
	remove_ext(c_tmp) ;
	remove_path(c_tmp) ;
	FILE *f,*f_tcl;
	sprintf(fout,"%s_VMD.sh",pdb_name);
	f = fopen(fout, "w") ;
	if(f){
		sprintf(fout2,"%s.tcl",pdb_name);

		f_tcl=fopen(fout2,"w");
		if(f_tcl){
			/* Write bash script for visualization using VMD */
			fprintf(f,"#!/bin/bash\nvmd %s -e %s.tcl\n",pdb_out_name,c_tmp);
                        fflush(f);
			fclose(f);
			
			/* Make tcl script executable, and Write tcl script */
			sprintf(sys_cmd,"chmod +x %s",fout);
			status = system(sys_cmd);

			
			fprintf(f_tcl,"proc highlighting { colorId representation id selection } {\n");
			fprintf(f_tcl,"   set id [[atomselect $id $selection] molid]\n");
			fprintf(f_tcl,"   puts \"highlighting $id\"\n");
   			fprintf(f_tcl,"   mol delrep 0 $id\n");
   			fprintf(f_tcl,"   mol representation $representation\n");
   			fprintf(f_tcl,"   mol color $colorId\n");
   			fprintf(f_tcl,"   mol selection $selection\n");
   			fprintf(f_tcl,"   mol addrep $id\n}\n\n");
			fprintf(f_tcl,"set repr \"Points 10\"\n");
			fprintf(f_tcl,"highlighting ResID \"Points 10\" 0 \"resname STP\"\n");
   			fprintf(f_tcl,"set id [[atomselect 0 \"protein\"] molid]\n");
   			fprintf(f_tcl,"puts \"highlighting $id\"\n");
   			fprintf(f_tcl,"mol representation \"Lines\"\n");
   			fprintf(f_tcl,"mol material \"Transparent\"\n");
   			fprintf(f_tcl,"mol color Element\n");
   			fprintf(f_tcl,"mol selection \"protein\"\n");
   			fprintf(f_tcl,"mol addrep $id\n");
   			fprintf(f_tcl,"set id [[atomselect 0 \"not protein and not resname STP\"] molid]\n");
   			fprintf(f_tcl,"puts \"highlighting $id\"\n");
   			fprintf(f_tcl,"mol representation \"Bonds\"\n");
   			fprintf(f_tcl,"mol color Element\n");
   			fprintf(f_tcl,"mol selection \"not protein and not resname STP\"\n");
   			fprintf(f_tcl,"mol addrep $id\n\n");
			
                        
                        fprintf(f_tcl,"mol new \"../%s.pdb\"\n",c_tmp);
                        fprintf(f_tcl,"mol selection \"not protein and not water\" \n \
                                mol material \"Opaque\" \n \
                                mol delrep 0 1 \n \
                                mol representation \"Lines 10\" \n \
                                mol addrep 1 \n \
                                highlighting Element \"NewCartoon\" 1 \"protein\"\n \
                                mol representation \"NewCartoon\" \n \
                                mol addrep $id \n \
                                mol new \"%s_pockets.pqr\"\n \
                                mol selection \"all\" \n \
                                mol material \"Glass1\" \n \
                                mol delrep 0 2 \n \
                                mol representation \"VDW\" \n \
                                mol color ResID 2 \n \
                                mol addrep 2 \n",c_tmp);
                        fclose(f_tcl);

		}
		else {
			fprintf(stderr, "! The file %s could not be opened!\n", fout2);
		}
		
	}
	else{
		fprintf(stderr, "! The file %s could not be opened!\n", fout);
	}

}

/**-----------------------------------------------------------------------------
   ## FUNCTION: 
	write_pymol
   -----------------------------------------------------------------------------
   ## SPECIFICATION: 
	write visualization script for PyMol
   -----------------------------------------------------------------------------
   ## PARAMETERS:
	@ char *pdb_name     : pdb code
	@ char *pdb_out_name : name of the pdb output file
   -----------------------------------------------------------------------------
   ## RETURN: void
   -----------------------------------------------------------------------------
*/
void write_pymol(char *pdb_name,char *pdb_out_name)
{
	char fout[250] = "" ;
	char fout2[250] = "" ;
	char sys_cmd[250] ="";
	FILE *f,*f_pml;
	char c_tmp[255];
	int status ;
	
	strcpy(c_tmp,pdb_name);
	remove_ext(c_tmp) ;
	remove_path(c_tmp) ;
	sprintf(fout,"%s_PYMOL.sh",pdb_name);
	f = fopen(fout, "w") ;
	if(f){
		sprintf(fout2,"%s.pml",pdb_name);
		f_pml=fopen(fout2,"w");
		if(f_pml){
			/* Write bash script for visualization using VMD */
			fprintf(f,"#!/bin/bash\npymol %s.pml\n",c_tmp);
                        fflush(f);
			fclose(f);
			
			sprintf(sys_cmd,"chmod +x %s",fout);
			status = system(sys_cmd);
			/* Write pml script */
			fprintf(f_pml,"from pymol import cmd,stored\n");
			fprintf(f_pml,"load %s\n",pdb_out_name);
			fprintf(f_pml,"#select pockets, resn STP\n");
   			fprintf(f_pml,"stored.list=[]\n");
   			fprintf(f_pml,"cmd.iterate(\"(resn STP)\",\"stored.list.append(resi)\")	#read info about residues STP\n");
   			fprintf(f_pml,"#print stored.list\n");
   			fprintf(f_pml,"lastSTP=stored.list[-1]	#get the index of the last residu\n");
   			fprintf(f_pml,"hide lines, resn STP\n\n");
			fprintf(f_pml,"#show spheres, resn STP\n");
			fprintf(f_pml,"for my_index in range(1,int(lastSTP)+1): cmd.select(\"pocket\"+str(my_index), \"resn STP and resi \"+str(my_index))\n");
			fprintf(f_pml,"for my_index in range(2,int(lastSTP)+2): cmd.color(my_index,\"pocket\"+str(my_index))\n");
   			fprintf(f_pml,"for my_index in range(1,int(lastSTP)+1): cmd.show(\"spheres\",\"pocket\"+str(my_index))\n");
   			fprintf(f_pml,"for my_index in range(1,int(lastSTP)+1): cmd.set(\"sphere_scale\",\"0.3\",\"pocket\"+str(my_index))\n");
   			fprintf(f_pml,"for my_index in range(1,int(lastSTP)+1): cmd.set(\"sphere_transparency\",\"0.1\",\"pocket\"+str(my_index))\n");
   			
			fclose(f_pml);
		}
		else {
			fprintf(stderr, "! The file %s could not be opened!\n", fout2);
		}
	}
	else{
		fprintf(stderr, "! The file %s could not be opened!\n", fout);
	}
}
