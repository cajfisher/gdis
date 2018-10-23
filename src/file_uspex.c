/*
Copyright (C) 2018 by Okadome Valencia

hubert.valencia _at_ imass.nagoya-u.ac.jp

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

The GNU GPL can also be found at http://www.gnu.org
*/

/* USPEX Parser:
 * supported versions VER 10.1, VER 9.4.4
 * All structure geometries are displayed "as an animation", additionally
 * also some graphs are produced depending on the USPEX calculation type:
 * Graph ALL    shows all structures energies vs generation number 
 * Graph BEST   shows the best _fitness_ structures per generation
 * Graph PATH   shows the image energies of VCNEB/TPS optimization
 * Graph COMP_X shows the structure energy for each X atomic ratio
 *
 * All graphs are selectable, ie. selecting a value represented by
 * a square on a graph will load the structure on the main display
 *
 * */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "error.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "zmatrix.h"
#include "zmatrix_pak.h"
#include "model.h"
#include "interface.h"
#include "graph.h"
#include "file_vasp.h"
#include "file_uspex.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

#define DEBUG_USPEX_READ 0

/***********************************/
/* Initialize uspex_calc structure */
/***********************************/
void init_uspex_parameters(uspex_calc_struct *uspex_calc){
	/*Initialize uspex_calc with model independant parameters.*/
	_UC.name=NULL;
	_UC.filename=NULL;
	_UC.path=NULL;
	_UC.special=0;
	_UC.calculationMethod=US_CM_USPEX;
	_UC.calculationType=300;
	_UC._calctype_dim=3;
	_UC._calctype_mol=FALSE;
	_UC._calctype_var=FALSE;
	_UC._calctype_mag=FALSE;	//VER 10.1
	_UC.optType=US_OT_ENTHALPY;
	_UC.new_optType=NULL;		//VER 10.1
	_UC.anti_opt=FALSE;
	_UC._nspecies=0;
	_UC.atomType=NULL;
	_UC._var_nspecies=0;
	_UC.numSpecies=NULL;
	_UC.magRatio[0]=0.1;		//VER 10.1
	_UC.magRatio[1]=0.225;		//VER 10.1
	_UC.magRatio[2]=0.225;		//VER 10.1
	_UC.magRatio[3]=0.225;		//VER 10.1
	_UC.magRatio[4]=0.225;		//VER 10.1
	_UC.magRatio[5]=0.0;		//VER 10.1
	_UC.magRatio[6]=0.0;		//VER 10.1
	_UC.ldaU=NULL;			//VER 10.1
	_UC.ExternalPressure=0.;
	_UC.valences=NULL;
	_UC.goodBonds=NULL;
	_UC.checkMolecules=TRUE;
	_UC.checkConnectivity=FALSE;
	_UC.fitLimit=0.;		//VER 10.1
	_UC.populationSize=0;
	_UC.initialPopSize=0;
	_UC.numGenerations=100;
	_UC.stopCrit=0;
	_UC.bestFrac=0.7;
	_UC.keepBestHM=0;
	_UC.reoptOld=FALSE;
	_UC.symmetries=NULL;
	_UC.fracGene=0.5;
	_UC.fracRand=0.2;
	_UC.fracTopRand=0.2;		//VER 10.1
	_UC.fracPerm=0.;
	_UC.fracAtomsMut=0.1;
	_UC.fracRotMut=0.;
	_UC.fracLatMut=0.;
	_UC.fracSpinMut=0.1;		//VER 10.1
	_UC.howManySwaps=0;
	_UC.specificSwaps=NULL;
	_UC.mutationDegree=0;
	_UC.mutationRate=0.5;
	_UC.DisplaceInLatmutation=1.0;
	_UC.AutoFrac=FALSE;
	_UC.minVectorLength=0.;
	_UC.IonDistances=NULL;
	_UC.constraint_enhancement=TRUE;
	_UC.MolCenters=NULL;
	_UC.Latticevalues=NULL;
	_UC.splitInto=NULL;
	_UC._num_opt_steps=0;
	_UC.abinitioCode=NULL;
	_UC.KresolStart=NULL;
	_UC.vacuumSize=NULL;
	_UC.numParallelCalcs=1;
	_UC.commandExecutable=NULL;
	_UC.whichCluster=0;
	_UC.remoteFolder=NULL;
	_UC.PhaseDiagram=FALSE;
	_UC.RmaxFing=10.0;
	_UC.deltaFing=0.08;
	_UC.sigmaFing=0.03;
	_UC.antiSeedsActivation=5000;
	_UC.antiSeedsMax=0.;
	_UC.antiSeedsSigma=0.001;
	_UC.doSpaceGroup=TRUE;
	_UC.SymTolerance=0.1;
	_UC.repeatForStatistics=1;
	_UC.stopFitness=0.;
	_UC.fixRndSeed=0;		//VER 10.1
	_UC.collectForces=FALSE;
	_UC.ordering_active=TRUE;
	_UC.symmetrize=FALSE;
	_UC.valenceElectr=NULL;
	_UC.percSliceShift=1.0;
	_UC.dynamicalBestHM=2;
	_UC.softMutOnly=NULL;
	_UC.maxDistHeredity=0.5;
	_UC.manyParents=0;
	_UC.minSlice=0.;
	_UC.maxSlice=0.;
	_UC.numberparents=2;
	_UC.BoltzTraP_T_max=800.0;	//VER 10.1
	_UC.BoltzTraP_T_delta=50.0;	//VER 10.1
	_UC.BoltzTraP_T_efcut=0.15;	//VER 10.1
	_UC.TE_T_interest=300.0;	//VER 10.1
	_UC.TE_threshold=0.5;		//VER 10.1
	_UC.TE_goal=US_BT_ZT;		//VER 10.1
	_UC.thicknessS=2.0;
	_UC.thicknessB=3.0;
	_UC.reconstruct=1;
	_UC.firstGeneMax=11;
	_UC.minAt=0;
	_UC.maxAt=0;
	_UC.fracTrans=0.1;
	_UC.howManyTrans=0.2;
	_UC.specificTrans=NULL;
	_UC.GaussianWidth=0.;
	_UC.GaussianHeight=0.;
	_UC.FullRelax=2;
	_UC.maxVectorLength=0.;
	_UC.PSO_softMut=1;
	_UC.PSO_BestStruc=1;
	_UC.PSO_BestEver=1;
	_UC.vcnebType=110;
	_UC._vcnebtype_method=1;
	_UC._vcnebtype_img_num=TRUE;
	_UC._vcnebtype_spring=FALSE;
	_UC.numImages=9;
	_UC.numSteps=200;
	_UC.optReadImages=2;
	_UC.optimizerType=1;
	_UC.optRelaxType=3;
	_UC.dt=0.05;
	_UC.ConvThreshold=0.003;
	_UC.VarPathLength=0.;
	_UC.K_min=5;
	_UC.K_max=5;
	_UC.Kconstant=5;
	_UC.optFreezing=FALSE;
	_UC.optMethodCIDI=0;
	_UC.startCIDIStep=100;
	_UC.pickupImages=NULL;
	_UC.FormatType=2;/*MUST BE VASP FORMAT!*/
	_UC.PrintStep=1;
	_UC.numIterations=1000;		//VER 10.1
	_UC.speciesSymbol=NULL;		//VER 10.1
	_UC.mass=NULL;			//VER 10.1
	_UC.amplitudeShoot[0]=0.1;	//VER 10.1
	_UC.amplitudeShoot[1]=0.1;	//VER 10.1
	_UC.magnitudeShoot[0]=1.05;	//VER 10.1
	_UC.magnitudeShoot[1]=1.05;	//VER 10.1
	_UC.shiftRatio=0.1;		//VER 10.1
	_UC.orderParaType=TRUE;		//VER 10.1 NOTE: here we select fingerprint method as default!
	_UC.opCriteria[0]=0.;		//VER 10.1
	_UC.opCriteria[1]=0.;		//VER 10.1
	_UC.cmdOrderParameter=NULL;			//VER 10.1
	_UC.cmdEnthalpyTemperature=NULL;		//VER 10.1
	_UC.orderParameterFile=g_strdup("fp.dat");	//VER 10.1
	_UC.enthalpyTemperatureFile=g_strdup("HT.dat");	//VER 10.1
	_UC.trajectoryFile=g_strdup("traj.dat");	//VER 10.1
	_UC.MDrestartFile=g_strdup("traj.restart");	//VER 10.1
}
void free_uspex_parameters(uspex_calc_struct *uspex_calc){
	/*free the sub-structures, and set it back to init*/
	g_free(_UC.name);
	g_free(_UC.path);
	g_free(_UC.filename);
	g_free(_UC.new_optType);	//VER 10.1
	g_free(_UC.atomType);
	g_free(_UC.numSpecies);
	g_free(_UC.ldaU);		//VER 10.1
	g_free(_UC.valences);
	g_free(_UC.goodBonds);
	g_free(_UC.symmetries);
	g_free(_UC.specificSwaps);
	g_free(_UC.IonDistances);
	g_free(_UC.MolCenters);
	g_free(_UC.Latticevalues);
	g_free(_UC.splitInto);
	g_free(_UC.abinitioCode);
	g_free(_UC.KresolStart);
	g_free(_UC.vacuumSize);
	g_free(_UC.commandExecutable);
	g_free(_UC.remoteFolder);
	g_free(_UC.valenceElectr);
	g_free(_UC.softMutOnly);
	g_free(_UC.specificTrans);
	g_free(_UC.pickupImages);
	g_free(_UC.speciesSymbol);		//VER 10.1
	g_free(_UC.mass);			//VER 10.1
	g_free(_UC.cmdOrderParameter);		//VER 10.1
	g_free(_UC.cmdEnthalpyTemperature);	//VER 10.1
	g_free(_UC.orderParameterFile);		//VER 10.1
	g_free(_UC.enthalpyTemperatureFile);	//VER 10.1
	g_free(_UC.trajectoryFile);		//VER 10.1
	g_free(_UC.MDrestartFile);		//VER 10.1
	init_uspex_parameters(uspex_calc);
}
/************************************************************************/
/* Read Parameters.txt and fill keywords of uspex_calc_struct structure */
/************************************************************************/
uspex_calc_struct *read_uspex_parameters(gchar *filename,gint safe_nspecies){
	FILE *vf;
	long int vfpos;
	gchar *line=NULL;
	gchar *ptr,*ptr2;
	gchar c;
	gint i,j,k;
	uspex_calc_struct *uspex_calc;
	/*tests*/
	if(filename==NULL) return NULL;
	vf = fopen(filename, "rt");
	if (!vf) return NULL;
	/**/
	uspex_calc = g_malloc(sizeof(uspex_calc_struct));
	init_uspex_parameters(uspex_calc);
/*this should be promoted*/
#define __Q(a) #a
/*some lazy defines*/
#define __GET_BOOL(value) if (find_in_string(__Q(value),line)!=NULL){\
	k=0;sscanf(line,"%i%*s",&(k));\
	_UC.value=(k==1);\
	g_free(line);line = file_read_line(vf);\
	continue;\
}
#define __GET_INT(value) if (find_in_string(__Q(value),line) != NULL) {\
	sscanf(line,"%i%*s",&(_UC.value));\
	g_free(line);line = file_read_line(vf);\
	continue;\
}
#define __GET_DOUBLE(value) if (find_in_string(__Q(value),line) != NULL) {\
	sscanf(line,"%lf%*s",&(_UC.value));\
	g_free(line);line = file_read_line(vf);\
	continue;\
}
#define __GET_STRING(value) if (find_in_string(__Q(value),line) != NULL) {\
	g_free(line);line = file_read_line(vf);\
	for(i=0;i<strlen(line);i++) if(line[i]=='\n') line[i]='\0';\
	_UC.value=g_strdup(line);\
	g_free(line);line = file_read_line(vf);\
	g_free(line);line = file_read_line(vf);\
	continue;\
}
#define __GET_CHARS(value) if (find_in_string(__Q(value),line) != NULL) {\
	ptr = g_strdup_printf("%%s : %s",__Q(value));\
	sscanf(line,ptr,&(_UC.value));\
	g_free(ptr);g_free(line);line = file_read_line(vf);\
	continue;\
}
#define __COUNT_NUM(pointer,count) while(*pointer!='\0'){\
	while((!g_ascii_isdigit(*pointer))&&(*pointer!='\0')) pointer++;\
	if(g_ascii_isdigit(*pointer)) count++;\
	while(g_ascii_isgraph(*pointer)&&(*pointer!='\0')) pointer++;\
}
#define __COUNT_ALNUM(pointer,count) while(*pointer!='\0'){\
	while((!g_ascii_isalnum(*pointer))&&(*pointer!='\0')) pointer++;\
	if(g_ascii_isalnum(*pointer)) count++;\
	while(g_ascii_isgraph(*pointer)&&(*pointer!='\0')) pointer++;\
}
#define __SKIP_BLANK(pointer) while(!g_ascii_isgraph(*pointer)&&(*ptr!='\0')) ptr++
	line = file_read_line(vf);
/* +++ 1st PASS: get important numbers*/
	while(line){
		if (find_in_string("atomType",line) !=NULL){
			g_free(line);line = file_read_line(vf);/*go next line*/
			ptr=&(line[0]);
			/*NEW: more safe*/
			_UC._nspecies=0;__COUNT_ALNUM(ptr,_UC._nspecies);
			if(_UC._nspecies<1) _UC._nspecies=1;/*there is at least 1 species*/
			g_free(line);line = file_read_line(vf);/*this is the EndAtomType line*/
			g_free(line);line = file_read_line(vf);
			continue;
		}
		if(find_in_string("abinitioCode",line)!=NULL){
		/*FIXME: in order to determine _num_opt_steps, *
		 * the number of items in abinitioCode is read *
		 * Even though 1 or VASP is stated as default. *
		 * ie. abinitioCode CANNOT be omited.  --ovhpa */
			g_free(line);line = file_read_line(vf);
			ptr=&(line[0]);
			/*NEW: more safe*/
			i=0;__COUNT_NUM(ptr,i);
			if(i>_UC._num_opt_steps) _UC._num_opt_steps=i;
			g_free(line);line = file_read_line(vf);/*this is the ENDabinit line*/
			g_free(line);line = file_read_line(vf);
			continue;
		}
		/*no catch*/
		g_free(line);line = file_read_line(vf);
	}
	g_free(line);
	if(_UC._nspecies==0){/*this one is mandatory*/
		if(safe_nspecies==0){
			line = g_strdup_printf("ERROR: Can't read USPEX output: no atom information!\n");
	                gui_text_show(ERROR, line);
	                g_free(line);
			fclose(vf);
			return(NULL);
		}else{
			/*this should *not* be the way to obtain _nspecies, atomType is mandatory!*/
			_UC._nspecies=safe_nspecies;
		}
	}
	rewind(vf);
/* +++ 2nd PASS: get everything else*/
	line = file_read_line(vf);
	while(line){
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: PROBE LINE: %s",line);
#endif
		if (find_in_string("calculationMethod",line) != NULL) {
			ptr=&(line[0]);
			__SKIP_BLANK(ptr);
			c=*ptr;//sscanf(line,"%c%*s",&(c));
			switch(c){
			case 'u':
			case 'U':/*USPEX method*/
				_UC.calculationMethod=US_CM_USPEX;
				break;
			case 'm':
			case 'M':
				/*MINHOP is unsupported*/
				ptr++;if((*ptr=='i')||(*ptr=='I')) _UC.calculationMethod=US_CM_UNKNOWN;
				else _UC.calculationMethod=US_CM_META;
				break;		//VER 10.1
			case 'v':
			case 'V':
				_UC.calculationMethod=US_CM_VCNEB;
				break;
			case 'p':
			case 'P':
				_UC.calculationMethod=US_CM_PSO;
				break;
			case 'T':
			case 't':
				_UC.calculationMethod=US_CM_TPS;
				break;		//VER 10.1
			case 'C':
			case 'c'://COPEX is unsupported VER 10.1
			default:
				_UC.calculationMethod=US_CM_UNKNOWN;
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("calculationType",line) != NULL) {
			ptr=&(line[0]);
			_UC.calculationType=0;
			__SKIP_BLANK(ptr);
			if(*ptr=='-') {/*2D-crystal*/
				j=-1;
				ptr++;
			}else j=1;
			if((*ptr=='s')||(*ptr=='S')) {/*VER 10.1: magnetic*/
				_UC._calctype_mag=TRUE;
				_UC.calculationType=1000;
				ptr++;
			}
			i=g_ascii_digit_value(*ptr);/*dimension*/
			if((i!=-2)&&((i<0)||(i>3))) {
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			_UC._calctype_dim=i;
			ptr++;
			i=g_ascii_digit_value(*ptr);/*molecular*/
			if((i<0)||(i>1)){
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			_UC._calctype_mol=(i==1);
			ptr++;
			i=g_ascii_digit_value(*ptr);/*variable*/
			if((i<0)||(i>1)){
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			_UC._calctype_var=(i==1);
			_UC.calculationType+=(_UC._calctype_dim*100)+(_UC._calctype_mol*10)+_UC._calctype_var;
			_UC.calculationType*=j;
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("optType",line) != NULL) {
			/*VER 10.1 has a new format !*/
			ptr=&(line[0]);
			__SKIP_BLANK(ptr);
			if(g_ascii_isalpha(*ptr)){
				if((*ptr=='m')||(*ptr=='M')){
					/*could be Min max*/
					j=0;
					if((*(ptr+1)=='a')&&(*(ptr+2)=='x')) {
						/*can be mag_moment -> checked x*/
						j=1;
						ptr+=4;
					}
					if(*(ptr+1)=='i') {
						/*can be only Min*/
						j=-1;
						ptr+=4;
					}
				}
				c=*ptr;
				switch(c){
				case 'E':
				case 'e':
					_UC.optType=US_OT_ENTHALPY;
					if(j==1) _UC.anti_opt=TRUE;
					break;
				case 'v':
				case 'V':
					_UC.optType=US_OT_VOLUME;
					if(j==1) _UC.anti_opt=TRUE;
					break;
				case 'h':
				case 'H':
					_UC.optType=US_OT_HARDNESS;
					if(j==-1) _UC.anti_opt=TRUE;
					break;
				case 's':
				case 'S':
					_UC.optType=US_OT_ORDER;
					if(j==-1) _UC.anti_opt=TRUE;
					break;
				case 'a':
				case 'A':
					_UC.optType=US_OT_DISTANCE;
					if(j==-1) _UC.anti_opt=TRUE;
					break;
				case 'd':
				case 'D':
					if((*(ptr+1)=='e')||(*(ptr+1)=='E')){
						/*VER 10.1 "density" replaces "aver_dist"*/
						_UC.optType=US_OT_DISTANCE;
						if(j==-1) _UC.anti_opt=TRUE;
						break;
					}else if((*(ptr+1)=='i')||(*(ptr+1)=='I')){
						if((*(ptr+5)=='s')||(*(ptr+5)=='S')) _UC.optType=US_OT_DIELEC_S;
						else if((*(ptr+5)=='g')||(*(ptr+5)=='G')) _UC.optType=US_OT_DIELEC_GAP;
						else _UC.optType=US_OT_UNKNOWN;
						if(j==-1) _UC.anti_opt=TRUE;
					} else _UC.optType=US_OT_UNKNOWN;
					break;
				case 'g':
				case 'G':
					_UC.optType=US_OT_GAP;
					break;
				case 'b':
				case 'B':
					if((*(ptr+1)=='a')||(*(ptr+1)=='A')){
						/*VER 10.1 "bandgap" replaces "gap"*/
						_UC.optType=US_OT_GAP;
						break;
					} else if((*(ptr+1)=='i')||(*(ptr+1)=='I')){
						_UC.optType=US_OT_2R;/*VER 10.1*/
						break;
					}
				case 'm':
				case 'M':
					_UC.optType=US_OT_MAG;
					break;
				case 'q':
				case 'Q':
					_UC.optType=US_OT_QE;
					break;
				case 'z':
				case 'Z':
					_UC.optType=US_OT_ZT;
					break;
				case 'f':
				case 'F':
					_UC.optType=US_OT_Fphon;
					break;
				default:
					_UC.optType=US_OT_UNKNOWN;
				}
			} else if(g_ascii_isdigit(*ptr)||(*ptr=='-')){
				/*load optType*/
				k=0;sscanf(ptr,"%i%*s",&(k));
				if(k<0) {
					_UC.anti_opt=TRUE;
					k*=-1;
				}
				switch (k){
				case 1:
					_UC.optType=US_OT_ENTHALPY;
					break;
				case 2:
					_UC.optType=US_OT_VOLUME;
					break;
				case 3:
					_UC.optType=US_OT_HARDNESS;
					break;
				case 4:
					_UC.optType=US_OT_ORDER;
					break;
				case 5:
					_UC.optType=US_OT_DISTANCE;
					break;
				case 6:
					_UC.optType=US_OT_DIELEC_S;
					break;
				case 7:
					_UC.optType=US_OT_GAP;
					break;
				case 8:
					_UC.optType=US_OT_DIELEC_GAP;
					break;
				case 9:
					_UC.optType=US_OT_MAG;
					break;
				case 10:
					_UC.optType=US_OT_QE;
					break;
				case 11:/*VER 10.1*/
					_UC.optType=US_OT_2R;
					break;
				case 14:/*VER 10.1*/
					_UC.optType=US_OT_ZT;
					break;
				case 17:/*VER 10.1*/
					_UC.optType=US_OT_Fphon;
					break;
				case 1101:
					_UC.optType=US_OT_BULK_M;
					break;
				case 1102:
					_UC.optType=US_OT_SHEAR_M;
					break;
				case 1103:
					_UC.optType=US_OT_YOUNG_M;
					break;
				case 1104:
					_UC.optType=US_OT_POISSON;
					break;
				case 1105:
					_UC.optType=US_OT_PUGH_R;
					break;
				case 1106:
					_UC.optType=US_OT_VICKERS_H;
					break;
				case 1107:
					_UC.optType=US_OT_FRACTURE;
					break;
				case 1108:
					_UC.optType=US_OT_DEBYE_T;
					break;
				case 1109:
					_UC.optType=US_OT_SOUND_V;
					break;
				case 1110:
					_UC.optType=US_OT_SWAVE_V;
					break;
				case 1111:
					_UC.optType=US_OT_PWAVE_V;
					break;
				default:
					_UC.optType=US_OT_UNKNOWN;
				}
			} else _UC.optType=US_OT_UNKNOWN;
			ptr=&(line[0]);
			i=0;__COUNT_ALNUM(ptr,i);
			if(i>1) {/*multiobjective optimization VER 10.1*/
				if(_UC.new_optType!=NULL) g_free(_UC.new_optType);
				_UC.new_optType=g_strdup(line);
			}else{
				if(_UC.new_optType!=NULL) g_free(_UC.new_optType);
				_UC.new_optType=NULL;
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("atomType",line) !=NULL){
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*look for atomType*/
			_UC.atomType = g_malloc(_UC._nspecies*sizeof(gint));
			for(i=0;i<_UC._nspecies;i++) _UC.atomType[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				__SKIP_BLANK(ptr);
				if(g_ascii_isdigit(*ptr)) _UC.atomType[i]=(gint)g_ascii_strtoull (ptr,NULL,10);/*user provided Z number*/
				else _UC.atomType[i]=elem_symbol_test(ptr);/*try user provided symbol*/
				if((_UC.atomType[i]<=0)||(_UC.atomType[i]>MAX_ELEMENTS-1)){/*invalid Z*/
					_UC.atomType[i]=0;/*allow it for now*/
				}
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndAtomType line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("numSpecies",line) != NULL) {
			vfpos=ftell(vf);/* flag */
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*1st get the total number of lines*/
			_UC._var_nspecies=0;
			do{
				g_free(line);line = file_read_line(vf);/*go next line*/
				_UC._var_nspecies++;
			}while(find_in_string("EndNumSpecies",line) == NULL);
			fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
			g_free(line);/*FIX: _VALGRIND_BUG_*/
			line = file_read_line(vf);/*first line of numSpecies*/
			_UC.numSpecies = g_malloc(_UC._nspecies*_UC._var_nspecies*sizeof(gint));
			for(i=0;i<_UC._nspecies*_UC._var_nspecies;i++) _UC.numSpecies[i]=0;
			for(j=0;j<_UC._var_nspecies;j++){
				ptr=&(line[0]);i=0;
				while((*ptr!='\n')&&(*ptr!='\0')){
					__SKIP_BLANK(ptr);
					_UC.numSpecies[i+j*_UC._nspecies]=(gint)g_ascii_strtoull(ptr,&ptr2,10);
					ptr=ptr2+1;
					i++;
				}
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("magRatio",line) != NULL) {
			/*a 7 line ratio for magnetic orders VER 10.1*/
			g_free(line);line = file_read_line(vf);/*go next line*/
			ptr=&(line[0]);
			i=0;ptr2=ptr;
			while((i<7)&&(*ptr!='\n')&&(*ptr!='\0')){
				/*two fromat are accepted: number or number/number*/
				__SKIP_BLANK(ptr);
				_UC.magRatio[i]=g_ascii_strtod(ptr,&ptr2);
				if(*ptr2=='/'){
					ptr=ptr2;
					_UC.magRatio[i]/=g_ascii_strtod(ptr,&ptr2);
				}
				ptr=ptr2+1;
				i++;
			}
			for(;i<7;i++) _UC.magRatio[i]=0.;/*fill missing values, if any*/
			g_free(line);
			line = file_read_line(vf);/*this is the EndMagRatio line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("ldaU",line) != NULL) {
			/*an array of U values for LDA+U VER 10.1*/
			g_free(line);line = file_read_line(vf);/*go next line*/
			ptr=&(line[0]);
			i=0;__COUNT_NUM(ptr,i);
			if(i>0){
				ptr=&(line[0]);ptr2=ptr;j=0;
				if(_UC.ldaU!=NULL) g_free(_UC.ldaU);
				_UC.ldaU=g_malloc(i*sizeof(gdouble));
				while((j<i)&&(*ptr!='\n')&&(*ptr!='\0')){
					__SKIP_BLANK(ptr);
					_UC.ldaU[j]=g_ascii_strtod(ptr,&ptr2);
					ptr=ptr2+1;
					j++;
				}
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndLdaU*/
			g_free(line);
			line = file_read_line(vf);
		}
		__GET_DOUBLE(ExternalPressure)
		if (find_in_string("valences",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*look for valences*/
			_UC.valences = g_malloc(_UC._nspecies*sizeof(gint));
			for(i=0;i<_UC._nspecies;i++) _UC.valences[i]=0;
			ptr=&(line[0]);ptr2=ptr;i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				__SKIP_BLANK(ptr);
				_UC.valences[i]=(gint)g_ascii_strtoull(ptr,&ptr2,10);
				ptr=ptr2+1;
				i++;
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndValences line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("goodBonds",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*look for goodBonds*/
			_UC.goodBonds = g_malloc(_UC._nspecies*_UC._nspecies*sizeof(gdouble));
			for(i=0;i<_UC._nspecies*_UC._nspecies;i++) _UC.goodBonds[i]=0.;
			j=0;
			while((j<_UC._nspecies)&&(find_in_string("EndGoodBonds",line)==NULL)){
				ptr=&(line[0]);ptr2=ptr;i=0;
				while((*ptr!='\n')&&(*ptr!='\0')){
					__SKIP_BLANK(ptr);
					_UC.goodBonds[i+j*_UC._nspecies]=g_ascii_strtod(ptr,&ptr2);
					ptr=ptr2+1;
					i++;
				}
				j++;
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_BOOL(checkMolecules);
		__GET_BOOL(checkConnectivity);
		__GET_DOUBLE(fitLimit);/*VER 10.1*/
		__GET_INT(populationSize);
		__GET_INT(initialPopSize);
		__GET_INT(numGenerations);
		__GET_INT(stopCrit);
		__GET_DOUBLE(bestFrac);
		__GET_INT(keepBestHM);
		__GET_BOOL(reoptOld);
		__GET_STRING(symmetries);
		__GET_DOUBLE(fracGene);
		__GET_DOUBLE(fracRand);
		__GET_DOUBLE(fracTopRand);/*VER 10.1*/
		__GET_DOUBLE(fracPerm);
		__GET_DOUBLE(fracAtomsMut);
		__GET_DOUBLE(fracRotMut);
		__GET_DOUBLE(fracLatMut);
		__GET_DOUBLE(fracSpinMut);/*VER 10.1*/
		__GET_INT(howManySwaps);
		__GET_STRING(specificSwaps);
		__GET_DOUBLE(mutationDegree);
		__GET_DOUBLE(mutationRate);
		__GET_DOUBLE(DisplaceInLatmutation);
		__GET_BOOL(AutoFrac);
		__GET_BOOL(minVectorLength);
		if (find_in_string("IonDistances",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*look for IonDistances*/
			_UC.IonDistances = g_malloc(_UC._nspecies*_UC._nspecies*sizeof(gdouble));
			for(i=0;i<_UC._nspecies*_UC._nspecies;i++) _UC.IonDistances[i]=0.;
			j=0;
			while((j<_UC._nspecies)&&(find_in_string("EndDistances",line)==NULL)){
				ptr=&(line[0]);ptr2=ptr;i=0;
				while((*ptr!='\n')&&(*ptr!='\0')){
					__SKIP_BLANK(ptr);
					_UC.IonDistances[i+j*_UC._nspecies]=g_ascii_strtod(ptr,&ptr2);
					ptr=ptr2+1;
					i++;
				}
				j++;
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(constraint_enhancement);
		if (find_in_string("MolCenters",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count number of molecules*/
			ptr=&(line[0]);
			_UC._nmolecules=0;__COUNT_ALNUM(ptr,_UC._nmolecules);
			if(_UC._nmolecules==0){
				_UC._nmolecules=1;/*there is at least one molecule*/
				_UC.MolCenters=g_malloc(1*sizeof(gdouble));
				_UC.MolCenters[0]=0.;/*this is bound to fail anyway*/
				g_free(line);
				line = file_read_line(vf);/*this is the EndMol line*/
			}else{
				_UC.MolCenters = g_malloc(_UC._nmolecules*_UC._nmolecules*sizeof(gdouble));
				for(i=0;i<_UC._nmolecules*_UC._nmolecules;i++) _UC.MolCenters[i]=0.;
				j=0;
				while((j<_UC._nmolecules)&&(find_in_string("EndMol",line)==NULL)){
					ptr=&(line[0]);ptr2=ptr;i=0;
					while((*ptr!='\n')&&(*ptr!='\0')){
						__SKIP_BLANK(ptr);
						_UC.MolCenters[i+j*_UC._nmolecules]=g_ascii_strtod(ptr,&ptr2);
						ptr=ptr2+1;
						i++;
					}
					j++;
					g_free(line);
					line = file_read_line(vf);
				}
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("Latticevalues",line) != NULL) {
			/*Latticevalues can be one or several volumes, one line crystal definition, or 3 lines lattice cell vectors...*/
			vfpos=ftell(vf);/* flag */
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count first line number of items*/
			ptr=&(line[0]);
			while(*ptr==' ') ptr++;
			_UC._nlatticevalues=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					_UC._nlatticevalues++;ptr++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}
			/*depending on varcomp*/
			if(_UC._calctype_var){
				/*varcomp calculation: should be only one line*/
				g_free(line);line = file_read_line(vf);/*go next line*/
				if(find_in_string("Endvalues",line) != NULL){/*only volumes values*/
					_UC.Latticevalues = g_malloc(_UC._nlatticevalues*sizeof(gdouble));
					for(i=0;i<_UC._nlatticevalues;i++) _UC.Latticevalues[i]=0.;
					fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
					g_free(line);line = file_read_line(vf);/*go next line*/
					ptr=&(line[0]);ptr2=ptr;i=0;
					while((*ptr!='\n')&&(*ptr!='\0')){
						__SKIP_BLANK(ptr);
						_UC.Latticevalues[i]=g_ascii_strtod(ptr,&ptr2);
						ptr=ptr2+1;
						i++;
					}
				}else{/*varcomp but more than one line -> bad setting*/
					g_free(line);
					line = g_strdup_printf("ERROR: USPEX bad Latticevalues volumes settings in varcomp Parameters.txt!\n");
					gui_text_show(ERROR, line);
					g_free(line);line = file_read_line(vf);/*go next line*/
				}
			}else{
				/*NOT varcomp calculation: should be a crystal definition*/
				if(_UC._nlatticevalues==6){/*crystal definition -- should we test for Endvalues?*/
					_UC.Latticevalues = g_malloc(_UC._nlatticevalues*sizeof(gdouble));
					for(i=0;i<_UC._nlatticevalues;i++) _UC.Latticevalues[i]=0.;
					fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
					g_free(line);line = file_read_line(vf);/*go next line*/
					ptr=&(line[0]);ptr2=ptr;i=0;
					while((*ptr!='\n')&&(*ptr!='\0')){
						__SKIP_BLANK(ptr);
						_UC.Latticevalues[i]=g_ascii_strtod(ptr,&ptr2);
						ptr=ptr2+1;
						i++;
					}
				}else if((_UC._nlatticevalues<=3)&&(_UC._nlatticevalues>1)){/*lattice vectors*/
					_UC.Latticevalues = g_malloc(_UC._nlatticevalues*_UC._nlatticevalues*sizeof(gdouble));
					for(i=0;i<_UC._nlatticevalues*_UC._nlatticevalues;i++) _UC.Latticevalues[i]=0.;
					fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
					g_free(line);line = file_read_line(vf);/*go next line*/
					ptr=&(line[0]);
					j=0;
					while((j<_UC._nlatticevalues)&&(find_in_string("Endvalues",line)==NULL)){
						ptr=&(line[0]);ptr2=ptr;i=0;
						while((*ptr!='\n')&&(*ptr!='\0')){
							__SKIP_BLANK(ptr);
							_UC.Latticevalues[i+j*_UC._nlatticevalues]=g_ascii_strtod(ptr,&ptr2);
							ptr=ptr2+1;
							i++;
						}
						j++;
						g_free(line);
						line = file_read_line(vf);
					}
				}else if(_UC._nlatticevalues==1){/*only the volume value*/
					_UC.Latticevalues = g_malloc(sizeof(gdouble));
					_UC.Latticevalues[0]=0.;
					_UC.Latticevalues[0]=g_ascii_strtod(ptr,NULL);
					g_free(line);
					line = file_read_line(vf);
				}else{/*bad setting*/
					g_free(line);
					line = g_strdup_printf("ERROR: USPEX bad Latticevalues setting in Parameters.txt!\n");
					gui_text_show(ERROR, line);
					g_free(line);line = file_read_line(vf);/*go next line*/
				}
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		/*TODO: update III*/
		if (find_in_string("splitInto",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count number of splits*/
			ptr=&(line[0]);
			while(*ptr==' ') ptr++;
			_UC._nsplits=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					_UC._nsplits++;ptr++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}
			/*look for splits*/
			_UC.splitInto = g_malloc(_UC._nsplits*sizeof(gint));
			for(i=0;i<_UC._nsplits;i++) _UC.splitInto[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%i%*s",&(_UC.splitInto[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndSplitInto line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("abinitioCode",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*look for abinitioCode*/
			_UC.abinitioCode = g_malloc(_UC._num_opt_steps*sizeof(gint));
			for(i=0;i<_UC._num_opt_steps;i++) _UC.abinitioCode[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if((*ptr==' ')||(*ptr=='(')) while((*ptr==' ')||(*ptr=='(')) ptr++;/*skip space(s) & parenthesis*/
				sscanf(ptr,"%i%*s",&(_UC.abinitioCode[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end (skip ')' if any)*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the ENDabinit line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("KresolStart",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*in my understanding, there is either 1 or _num_opt_steps values of KresolStart*/
			_UC.KresolStart = g_malloc(_UC._num_opt_steps*sizeof(gdouble));
			for(i=0;i<_UC._num_opt_steps;i++) _UC.KresolStart[i]=0.;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%lf%*s",&(_UC.KresolStart[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the Kresolend line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		if (find_in_string("vacuumSize",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*in my understanding, there is BETWEEN 1 to _num_opt_steps values of vacuumSize*/
			_UC.vacuumSize = g_malloc(_UC._num_opt_steps*sizeof(gdouble));
			for(i=0;i<_UC._num_opt_steps;i++) _UC.vacuumSize[i]=10.;/*implicit is 10.*/
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%lf%*s",&(_UC.vacuumSize[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the endVacuumSize line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(numParallelCalcs);
		if (find_in_string("commandExecutable",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*there is also 1<x<_num_opt_steps lines of commandExecutable?*/
			j=1;
			ptr = g_strdup(line);
			_UC.commandExecutable = ptr;
			g_free(line);
			line = file_read_line(vf);
			while((j<_UC._num_opt_steps)&&(find_in_string("EndExecutable",line)==NULL)){
				_UC.commandExecutable = g_strdup_printf("%s\n%s",ptr,line);
				g_free(ptr);
				g_free(line);
				line = file_read_line(vf);
				ptr = _UC.commandExecutable;
				j++;
			}
			/*remove last '\n'*/
			j = strlen(_UC.commandExecutable);
			ptr = &(_UC.commandExecutable[j]);
			while((*ptr!='\n')&&(j>0)) j--;
			ptr[j-1]='\0';
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(whichCluster);
		__GET_CHARS(remoteFolder);
		__GET_BOOL(PhaseDiagram);
		__GET_DOUBLE(RmaxFing);
		__GET_DOUBLE(deltaFing);
		__GET_DOUBLE(sigmaFing);
		__GET_INT(antiSeedsActivation);
		__GET_DOUBLE(antiSeedsMax);
		__GET_DOUBLE(antiSeedsSigma);
		__GET_BOOL(doSpaceGroup);
		if (find_in_string("SymTolerance",line) != NULL) {
			/*can be either a number or quantifier*/
			if(g_ascii_isalpha(line[0])){
				/*Either high, medium, or low*/
				switch(line[0]){
				case 'H':
				case 'h':
					_UC.SymTolerance=0.05;
					break;
				case 'M':
				case 'm':
					_UC.SymTolerance=0.10;
					break;
				case 'L':
				case 'l':
					_UC.SymTolerance=0.20;
					break;
				default:/*everything else is a bad setting*/
					g_free(line);
					line = g_strdup_printf("ERROR: USPEX bad SymTolerance setting in Parameters.txt!\n");
					gui_text_show(ERROR, line);
				}
			}else{
				/*get the number directly*/
				sscanf(line,"%lf%*s",&(_UC.SymTolerance));
			}
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(repeatForStatistics);
		__GET_DOUBLE(stopFitness);
		__GET_INT(fixRndSeed);/*VER 10.1*/
		__GET_BOOL(collectForces);
		__GET_BOOL(ordering_active);
		__GET_BOOL(symmetrize);
		if (find_in_string("valenceElectr",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*look for valenceElectr*/
			_UC.valenceElectr = g_malloc(_UC._nspecies*sizeof(gint));
			for(i=0;i<_UC._nspecies;i++) _UC.valenceElectr[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%i%*s",&(_UC.valenceElectr[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the endValenceElectr line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_DOUBLE(percSliceShift);
		__GET_INT(dynamicalBestHM);
		__GET_STRING(softMutOnly);/*TODO: rewrite as "smart" int array*/
		__GET_DOUBLE(maxDistHeredity);
		__GET_INT(manyParents);
		__GET_DOUBLE(minSlice);
		__GET_DOUBLE(maxSlice);
		__GET_INT(numberparents);
		__GET_DOUBLE(thicknessS);
		__GET_DOUBLE(thicknessB);
		__GET_INT(reconstruct);
		__GET_INT(firstGeneMax);
		__GET_INT(minAt);
		__GET_INT(maxAt);
		__GET_DOUBLE(fracTrans);
		__GET_DOUBLE(howManyTrans);
		if (find_in_string("specificTrans",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count number of specificTrans*/
			ptr=&(line[0]);
			_UC._nspetrans=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					_UC._nspetrans++;ptr++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}
			/*look for specificTrans*/
			_UC.specificTrans = g_malloc(_UC._nspetrans*sizeof(gint));
			for(i=0;i<_UC._nspetrans;i++) _UC.specificTrans[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%i%*s",&(_UC.specificTrans[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndTransSpecific line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_DOUBLE(GaussianWidth);
		__GET_DOUBLE(GaussianHeight);
		__GET_INT(FullRelax);
		__GET_DOUBLE(maxVectorLength);
		__GET_DOUBLE(PSO_softMut);
		__GET_DOUBLE(PSO_BestStruc);
		__GET_DOUBLE(PSO_BestEver);
		if (find_in_string("vcnebType",line) != NULL) {
			k=0;sscanf(line,"%i%*s",&(k));
			j=k;
			i=((int)(j/100))%10;
			if((i<1)||(i>2)) {
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			_UC._vcnebtype_method=i;
			j-=i*100;
			i=((int)(j/10))%10;
			if((i<0)||(i>1)){
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			_UC._vcnebtype_img_num=(i==1);
			j-=i*10;
			i=j%10;
			if((i<0)||(i>1)){
				g_free(line);
				line = file_read_line(vf);
				continue;
			}
			_UC._vcnebtype_spring=(i==1);
			_UC.vcnebType=k;
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(numImages);
		__GET_INT(numSteps);
		__GET_INT(optReadImages);
		__GET_INT(optimizerType);
		__GET_INT(optRelaxType);
		__GET_DOUBLE(dt);
		__GET_DOUBLE(ConvThreshold);
		__GET_DOUBLE(VarPathLength);
		__GET_DOUBLE(K_min);
		__GET_DOUBLE(K_max);
		__GET_DOUBLE(Kconstant);
		__GET_BOOL(optFreezing);
		__GET_INT(optMethodCIDI);
		__GET_INT(startCIDIStep);
		if (find_in_string("pickupImages",line) != NULL) {
			g_free(line);line = file_read_line(vf);/*go next line*/
			/*count number of picked-up images*/
			ptr=&(line[0]);
			_UC._npickimg=1;
			while(*ptr!='\0'){
				if(*ptr==' ') {
					_UC._npickimg++;ptr++;
					while(g_ascii_isgraph(*ptr)) ptr++;
				}else ptr++;
			}
			/*look for pickupImages*/
			_UC.pickupImages = g_malloc(_UC._npickimg*sizeof(gint));
			for(i=0;i<_UC._npickimg;i++) _UC.pickupImages[i]=0;
			ptr=&(line[0]);i=0;
			while((*ptr!='\n')&&(*ptr!='\0')){
				if(*ptr==' ') while(*ptr==' ') ptr++;/*skip space(s)*/
				sscanf(ptr,"%i%*s",&(_UC.pickupImages[i]));
				ptr++;i++;
				while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
			}
			g_free(line);
			line = file_read_line(vf);/*this is the EndPickupImages line*/
			g_free(line);
			line = file_read_line(vf);
			continue;
		}
		__GET_INT(FormatType);
		__GET_INT(PrintStep);
		/*in case no match was done:*/
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: PROBE FAIL LINE: %s",line);
#endif
		g_free(line);
		line = file_read_line(vf);
	}
	fclose(vf);
	return uspex_calc;
#undef __GET_BOOL
#undef __GET_INT
#undef __GET_DOUBLE
#undef __GET_STRING
#undef __GET_CHARS
}
void copy_uspex_parameters(uspex_calc_struct *src,uspex_calc_struct *dest){
	gint i;
#define _SRC (*src)
#define _DEST (*dest)
#define _CP(value) _DEST.value=_SRC.value
#define _COPY(value,size,type) do{\
	if((_SRC.value)!=NULL){\
		if(_DEST.value!=NULL) g_free(_DEST.value);\
		_DEST.value = g_malloc(size*sizeof(type));\
		memcpy(_DEST.value,_SRC.value,size*sizeof(type));\
	}else{\
		_DEST.value=NULL;\
	}\
}while(0)
#define _DBLCP(value,size) do{\
	if((_SRC.value)!=NULL){\
		if(_DEST.value!=NULL) g_free(_DEST.value);\
		_DEST.value = g_malloc(size*sizeof(gdouble));\
		for(i=0;i<size;i++) _DEST.value[i]=_SRC.value[i];\
	}else{\
		_DEST.value=NULL;\
	}\
}while(0)
#define _STRCP(value) do{\
	if((_SRC.value)!=NULL){\
		if(_DEST.value!=NULL) g_free(_DEST.value);\
		_DEST.value = g_strdup(_SRC.value);\
	}else{\
		_DEST.value=NULL;\
	}\
}while(0)
	/*1st free / init the destination*/
	free_uspex_parameters(dest);
	init_uspex_parameters(dest);
	/*copy*/
	_CP(name);
	_STRCP(filename);
	_STRCP(path);
	_CP(special);
	_CP(calculationMethod);
	_CP(calculationType);
	_CP(_calctype_dim);
	_CP(_calctype_mol);
	_CP(_calctype_var);
	_CP(optType);
	_CP(anti_opt);
	_CP(_nspecies);
	_COPY(atomType,_SRC._nspecies,gint);
	_CP(_var_nspecies);
	_COPY(numSpecies,(_SRC._nspecies*_SRC._var_nspecies),gint);
	_CP(ExternalPressure);
	_COPY(valences,_SRC._nspecies,gint);
	_DBLCP(goodBonds,_SRC._nspecies*_SRC._nspecies);
	_CP(checkMolecules);
	_CP(checkConnectivity);
	_CP(populationSize);
	_CP(initialPopSize);
	_CP(numGenerations);
	_CP(stopCrit);
	_CP(bestFrac);
	_CP(keepBestHM);
	_CP(reoptOld);
	_STRCP(symmetries);
	_CP(fracGene);
	_CP(fracRand);
	_CP(fracPerm);
	_CP(fracAtomsMut);
	_CP(fracRotMut);
	_CP(fracLatMut);
	_CP(howManySwaps);
	_STRCP(specificSwaps);
	_CP(mutationDegree);
	_CP(mutationRate);
	_CP(DisplaceInLatmutation);
	_CP(AutoFrac);
	_CP(minVectorLength);
	_DBLCP(IonDistances,_SRC._nspecies*_SRC._nspecies);
	_CP(constraint_enhancement);
	_CP(_nmolecules);
	_DBLCP(MolCenters,_SRC._nmolecules);
	if(_SRC._calctype_var){
		_DBLCP(Latticevalues,_SRC._nlatticevalues);
	}else{
		if(_SRC._nlatticevalues==6) _DBLCP(Latticevalues,6);
		else if(_SRC._nlatticevalues<=3) _DBLCP(Latticevalues,_SRC._nlatticevalues*_SRC._nlatticevalues);
		else if(_SRC._nlatticevalues==1) _DBLCP(Latticevalues,1);
	}
	_COPY(splitInto,_SRC._nsplits,gint);
	_COPY(abinitioCode,_SRC._num_opt_steps,gint);
	_DBLCP(KresolStart,_SRC._num_opt_steps);
	_DBLCP(vacuumSize,_SRC._num_opt_steps);
	_CP(numParallelCalcs);
	_STRCP(commandExecutable);/*unsure*/
	_CP(whichCluster);
	_STRCP(remoteFolder);
	_CP(PhaseDiagram);
	_CP(RmaxFing);
	_CP(deltaFing);
	_CP(sigmaFing);
	_CP(antiSeedsActivation);
	_CP(antiSeedsMax);
	_CP(antiSeedsSigma);
	_CP(doSpaceGroup);
	_CP(SymTolerance);
	_CP(repeatForStatistics);
	_CP(stopFitness);
	_CP(collectForces);
	_CP(ordering_active);
	_CP(symmetrize);
	_COPY(valenceElectr,_SRC._nspecies,gint);
	_CP(percSliceShift);
	_CP(dynamicalBestHM);
	_STRCP(softMutOnly);
	_CP(maxDistHeredity);
	_CP(manyParents);
	_CP(minSlice);
	_CP(maxSlice);
	_CP(numberparents);
	_CP(thicknessS);
	_CP(thicknessB);
	_CP(reconstruct);
	_CP(firstGeneMax);
	_CP(minAt);
	_CP(maxAt);
	_CP(fracTrans);
	_CP(howManyTrans);
	_COPY(specificTrans,_SRC._nspetrans,gint);
	_CP(GaussianWidth);
	_CP(GaussianHeight);
	_CP(FullRelax);
	_CP(maxVectorLength);
	_CP(PSO_softMut);
	_CP(PSO_BestStruc);
	_CP(PSO_BestEver);
	_CP(vcnebType);
	_CP(_vcnebtype_method);
	_CP(_vcnebtype_img_num);
	_CP(_vcnebtype_spring);
	_CP(numImages);
	_CP(numSteps);
	_CP(optReadImages);
	_CP(optimizerType);
	_CP(optRelaxType);
	_CP(dt);
	_CP(ConvThreshold);
	_CP(VarPathLength);
	_CP(K_min);
	_CP(K_max);
	_CP(Kconstant);
	_CP(optFreezing);
	_CP(optMethodCIDI);
	_CP(startCIDIStep);
	_COPY(pickupImages,_SRC._npickimg,gint);
	_CP(FormatType);
	_CP(PrintStep);
}

/*****************************************/
/* OUTPUT a minimal uspex parameter file */
/*****************************************/
gint dump_uspex_parameters(gchar *filename,uspex_calc_struct *uspex_calc){
#define __TITLE(caption) do{\
	gchar *tmp=g_strnfill (strlen(caption),'*');\
	fprintf(vf,"%s\n",tmp);\
	fprintf(vf,"%s\n",caption);\
	fprintf(vf,"%s\n",tmp);\
	g_free(tmp);\
}while(0)
#define __OUT_BOOL(value) do{\
	if(_UC.value==1) fprintf(vf,"1\t: %s\n",__Q(value));\
	else fprintf(vf,"0\t: %s\n",__Q(value));\
	is_w++;\
}while(0)
#define __OUT_INT(value) do{\
	fprintf(vf,"%i\t: %s\n",_UC.value,__Q(value));\
	is_w++;\
}while(0)
#define __OUT_DOUBLE(value) do{\
	fprintf(vf,"%.5f\t: %s\n",_UC.value,__Q(value));\
	is_w++;\
}while(0)
#define __OUT_STRING(value) do{\
	fprintf(vf,"%s\t: %s\n",_UC.value,__Q(value));\
	is_w++;\
}while(0)
/*because there is no constant logic in end_tag, we need to supply it*/
#define __OUT_BK_INT(value,end_tag,number) do{\
	fprintf(vf,"%% %s\n",__Q(value));\
	for(i=0;i<(number);i++) fprintf(vf,"%i ",_UC.value[i]);\
	fprintf(vf,"\n");fprintf(vf,"%% %s\n",end_tag);\
	is_w++;\
}while(0)
#define __OUT_BK_DOUBLE(value,end_tag,number) do{\
        fprintf(vf,"%% %s\n",__Q(value));\
        for(i=0;i<(number);i++) fprintf(vf,"%.5f ",_UC.value[i]);\
        fprintf(vf,"\n");fprintf(vf,"%% %s\n",end_tag);\
	is_w++;\
}while(0)
#define __OUT_BK_STRING(value,end_tag) do{\
	fprintf(vf,"%% %s\n",__Q(value));\
	fprintf(vf,"%s\n",_UC.value);\
	fprintf(vf,"%% %s\n",end_tag);\
	is_w++;\
}while(0)
#define __OUT_TMAT_DOUBLE(value,end_tag,number) do{\
	fprintf(vf,"%% %s\n",__Q(value));\
	for(i=0;i<(number);i++) {\
		for(j=0;j<(number);j++) fprintf(vf,"%.5f ",_UC.value[j+i*(number)]);\
		fprintf(vf,"\n");\
	}\
	fprintf(vf,"%% %s\n",end_tag);\
	is_w++;\
}while(0)

	/**/
        FILE *vf,*dest;
        long int vfpos;
	gint is_w=0;
        gchar *line=NULL;
        gint i,j;//,k;
        /*tests*/
        if(filename==NULL) return -1;
        if(uspex_calc==NULL) return -1;
	/*check integrity of uspex_calc structure.
	 *ie: check if mandatory settings are set.*/
	if(_UC.atomType==NULL) {
		line = g_strdup_printf("ERROR: USPEX - missing atomType!\n");
		gui_text_show(ERROR, line);g_free(line);
		return -5;
	}
	if(_UC.numSpecies==NULL) {
		line = g_strdup_printf("ERROR: USPEX - missing numSpecies!\n");
		gui_text_show(ERROR, line);g_free(line);
		return -5;
	}
	if((_UC.ExternalPressure==0.)&&(_UC.calculationMethod==US_CM_META)) {
		line = g_strdup_printf("ERROR: USPEX - META calculation but ExternalPressure is 0.!\n");
		gui_text_show(ERROR, line);g_free(line);
		return -5;
	}
	if((_UC._calctype_mol)&&(_UC._nmolecules>1)&&(_UC.MolCenters==NULL)) {
		line = g_strdup_printf("ERROR: USPEX - molecular calculation but missing MolCenters!\n");
		gui_text_show(ERROR, line);g_free(line);
		return -5;
	}
	if((_UC.repeatForStatistics>1)&&(_UC.stopFitness==0.)){
		line = g_strdup_printf("ERROR: USPEX - repeatForStatistics is >1 but stopFitness is 0.!\n");
		gui_text_show(ERROR, line);g_free(line);
		return -5;
	}
	if((_UC._calctype_var)&&(_UC.minAt==0)){
		line = g_strdup_printf("ERROR: USPEX - varcomp calculation but missing minAt!\n");
		gui_text_show(ERROR, line);g_free(line);
		return -5;
	}
	if((_UC._calctype_var)&&(_UC.maxAt==0)){
		line = g_strdup_printf("ERROR: USPEX - varcomp calculation but missing maxAt!\n");
		gui_text_show(ERROR, line);g_free(line);
		return -5;
	}
	if((_UC.calculationMethod==US_CM_META)&&(_UC.maxVectorLength=0.)){
		line = g_strdup_printf("ERROR: USPEX - META calculation but missong maxVectorLength!\n");
		gui_text_show(ERROR, line);g_free(line);
		return -5;
	}
	/*open output*/
	dest = fopen(filename, "w");
	if (!dest) return -2;
	vf = tmpfile();
	if(!vf) {
		fclose(dest);
		return -3;
	}
        /*output resulting "Parameter.txt" input*/
	fprintf(vf,"PARAMETERS EVOLUTIONARY ALGORITHM\n");
	fprintf(vf,"GENERATED BY GDIS %4.2f.%d (C) %d\n",VERSION,PATCH,YEAR);
	if(_UC.name==NULL) _UC.name=g_strdup("(unnamed)");
	line=g_strdup_printf("* CALCULATION: %s *",_UC.name);
	__TITLE(line);
	g_free(line);
	line=g_strdup_printf("*       TYPE OF RUN AND SYSTEM       *");
	__TITLE(line);
	g_free(line);
	/*always print:*/
	switch(_UC.calculationMethod){
	case US_CM_USPEX:
		fprintf(vf,"USPEX\t: calculationMethod\n");
		break;
	case US_CM_META:
		fprintf(vf,"META\t: calculationMethod\n");
		break;
	case US_CM_VCNEB:
		fprintf(vf,"VCNEB\t: calculationMethod\n");
		break;
	case US_CM_PSO:
		fprintf(vf,"PSO\t: calculationMethod\n");
		break;
	case US_CM_UNKNOWN:
	default:
		fprintf(vf,"???\t: calculationMethod (unsupported method)\n");
	}
	__OUT_INT(calculationType);
	if(_UC.anti_opt) __OUT_INT(optType);
	else {
		fprintf(vf,"%i\t: optType\n",-1*_UC.optType);
	}
	__OUT_BK_INT(atomType,"EndAtomType",_UC._nspecies);
	if(_UC._var_nspecies==1){
		__OUT_BK_INT(numSpecies,"EndNumSpecies",_UC._nspecies);
		is_w++;
	}else{
		fprintf(vf,"%% numSpecies\n");
		for(i=0;i<_UC._var_nspecies;i++){
			for(j=0;j<_UC._nspecies;j++){
				fprintf(vf,"%i ",_UC.numSpecies[j+i*_UC._nspecies]);
			}
			fprintf(vf,"\n");
		}
		fprintf(vf,"%% EndNumSpecies\n");
		is_w++;
	}
	if(_UC.calculationMethod==US_CM_META) __OUT_DOUBLE(ExternalPressure);
	if(_UC.ExternalPressure!=0.) __OUT_DOUBLE(ExternalPressure);
	/*print when NOT default or unset*/
	if(_UC.valences!=NULL) __OUT_BK_INT(valences,"endValences",_UC._nspecies);
	if(_UC.goodBonds!=NULL) __OUT_TMAT_DOUBLE(goodBonds,"EndGoodBonds",_UC._nspecies);/*TODO: use of a single value is unsupported*/
	if(!_UC.checkMolecules) __OUT_BOOL(checkMolecules);
	if(_UC.checkConnectivity) __OUT_BOOL(checkConnectivity);
vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       POPULATION       *");
	__TITLE(line);
	g_free(line);
	if(_UC.populationSize!=0) __OUT_INT(populationSize);
	if((_UC.initialPopSize!=0)&&(_UC.initialPopSize!=_UC.populationSize)) __OUT_INT(initialPopSize);
	if(_UC.numGenerations!=100) __OUT_INT(numGenerations);
	if(_UC.stopCrit!=0) __OUT_INT(stopCrit);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       SURVIVAL OF THE FITTEST & SELECTION       *");
	__TITLE(line);
	g_free(line);
	if(_UC.bestFrac!=0.7) __OUT_DOUBLE(bestFrac);
	if(_UC.keepBestHM!=0) __OUT_INT(keepBestHM);
	if(_UC.reoptOld) __OUT_BOOL(reoptOld);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;	
	line=g_strdup_printf("*       STRUCTURE GENERATION & VARIATION OPERATORS       *");
	__TITLE(line);
	g_free(line);
	if(_UC.symmetries!=NULL) __OUT_BK_STRING(symmetries,"endSymmetries");
	if(_UC.fracGene!=0.5) __OUT_DOUBLE(fracGene);
	if(_UC.fracRand!=0.2) __OUT_DOUBLE(fracRand);
	if((_UC._nspecies>1)&&(_UC.fracPerm!=0.1)) __OUT_DOUBLE(fracPerm);
	if((_UC._nspecies==1)&&(_UC.fracPerm!=0.)) __OUT_DOUBLE(fracPerm);
	if(_UC.fracAtomsMut!=0.1) __OUT_DOUBLE(fracAtomsMut);
	if((_UC._calctype_dim==3)&&(_UC._calctype_mol)){
		if(_UC.fracRotMut!=0.1) __OUT_DOUBLE(fracRotMut);
	}else{
		if(_UC.fracRotMut!=0.) __OUT_DOUBLE(fracRotMut);
	}
	if((_UC.FullRelax==0)||(_UC.optRelaxType==1)){
		if(_UC.fracLatMut!=0.) __OUT_DOUBLE(fracLatMut);
	}else{
		if(_UC.fracLatMut!=0.1) __OUT_DOUBLE(fracLatMut);
	}
	if(_UC.howManySwaps!=0) __OUT_INT(howManySwaps);
	if(_UC.specificSwaps!=NULL) __OUT_BK_STRING(specificSwaps,"EndSpecific");
	if(_UC.mutationDegree!=0.) __OUT_DOUBLE(mutationDegree);
	if(_UC.mutationRate!=0.5) __OUT_DOUBLE(mutationRate);
	if(_UC.DisplaceInLatmutation!=1.0) __OUT_DOUBLE(DisplaceInLatmutation);
	if(_UC.AutoFrac) __OUT_BOOL(AutoFrac);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       CONSTRAINTS       *");
	__TITLE(line);
	g_free(line);
	if(_UC.minVectorLength!=0.) __OUT_DOUBLE(minVectorLength);
	if(_UC.IonDistances!=NULL) __OUT_TMAT_DOUBLE(IonDistances,"EndDistances",_UC._nspecies);
	if(!_UC.constraint_enhancement) __OUT_BOOL(constraint_enhancement);
	if(_UC.MolCenters!=NULL) __OUT_TMAT_DOUBLE(MolCenters,"EndMol",_UC._nmolecules);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       CELL       *");
	__TITLE(line);
	g_free(line);
	if(_UC.Latticevalues!=NULL){
		if(_UC._calctype_var) {
			/*should be 1 line*/
			__OUT_BK_DOUBLE(Latticevalues,"Endvalues",_UC._nlatticevalues);
		}else{/*more than one line if 1<_UC._nlatticevalues<=3*/
			if((_UC._nlatticevalues>1)&&(_UC._nlatticevalues<=3)) __OUT_TMAT_DOUBLE(Latticevalues,"Endvalues",_UC._nlatticevalues);
			else __OUT_BK_DOUBLE(Latticevalues,"Endvalues",_UC._nlatticevalues);
		}
	}
	if(_UC.splitInto!=NULL) __OUT_BK_INT(splitInto,"EndSplitInto",_UC._nsplits);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       AB INITIO CALCULATION       *");
	__TITLE(line);
	g_free(line);
	if(_UC.abinitioCode!=NULL) __OUT_BK_INT(abinitioCode,"ENDabinit",_UC._num_opt_steps);
	if(_UC.KresolStart!=NULL) __OUT_BK_DOUBLE(KresolStart,"Kresolend",_UC._num_opt_steps);
	if(_UC.vacuumSize!=NULL) __OUT_BK_DOUBLE(vacuumSize,"endVacuumSize",_UC._num_opt_steps);
	if(_UC.numParallelCalcs!=1) __OUT_INT(numParallelCalcs);
	__OUT_BK_STRING(commandExecutable,"EndExecutable");/*mandatory block*/
	if(_UC.whichCluster!=0) __OUT_INT(whichCluster);
	if((_UC.whichCluster==2)&&(_UC.remoteFolder!=NULL)) __OUT_STRING(remoteFolder);
	if(_UC.PhaseDiagram) __OUT_BOOL(PhaseDiagram);
/*note that mandatory block will always make is_w>0*/
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       FINGERPRINT SETTINGS       *");
	__TITLE(line);
	g_free(line);
	if(_UC.RmaxFing!=10.) __OUT_DOUBLE(RmaxFing);
	if(_UC.deltaFing!=0.08) __OUT_DOUBLE(deltaFing);
	if(_UC.sigmaFing!=0.03) __OUT_DOUBLE(sigmaFing);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       ANTISEEDS SETTINGS       *");
	__TITLE(line);
	g_free(line);
	if(_UC.antiSeedsActivation!=5000) __OUT_INT(antiSeedsActivation);
	if(_UC.antiSeedsMax!=0.) __OUT_DOUBLE(antiSeedsMax);
	if(_UC.antiSeedsSigma!=0.001) __OUT_DOUBLE(antiSeedsSigma);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       SPACE GROUP DETERMINATION       *");
	__TITLE(line);
	g_free(line);
	if((_UC._calctype_dim==3)&&(!_UC.doSpaceGroup)) __OUT_BOOL(doSpaceGroup);
	if((_UC._calctype_dim!=3)&&(_UC.doSpaceGroup)) __OUT_BOOL(doSpaceGroup);
	if(_UC.SymTolerance!=0.1) __OUT_DOUBLE(SymTolerance);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       DEVELOPERS       *");	
	__TITLE(line);
	g_free(line);
	if(_UC.repeatForStatistics!=1) {
		__OUT_INT(repeatForStatistics);
		if(_UC.stopFitness!=0.) __OUT_DOUBLE(stopFitness);
	}
	if(_UC.collectForces) __OUT_BOOL(collectForces);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       SELDOM       *");
	__TITLE(line);
	g_free(line);
	if(!_UC.ordering_active) __OUT_BOOL(ordering_active);
	if(_UC.symmetrize) __OUT_BOOL(symmetrize);
	if(_UC.valenceElectr!=NULL) __OUT_BK_INT(valenceElectr,"endValenceElectr",_UC._nspecies);
	if(_UC.percSliceShift!=1.0) __OUT_DOUBLE(percSliceShift);
	if(_UC.dynamicalBestHM!=2) __OUT_INT(dynamicalBestHM);
	if(_UC.softMutOnly!=NULL) __OUT_BK_STRING(softMutOnly,"EndSoftOnly");
	if(_UC.maxDistHeredity!=0.5) __OUT_DOUBLE(maxDistHeredity);
	if(_UC.manyParents!=0) __OUT_INT(manyParents);
	if(_UC.minSlice!=0.) __OUT_DOUBLE(minSlice);
	if(_UC.maxSlice!=0.) __OUT_DOUBLE(maxSlice);
	if(_UC.numberparents!=2) __OUT_INT(numberparents);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       SURFACES       *");
	__TITLE(line);
	g_free(line);
	if(_UC.thicknessS!=2.0) __OUT_DOUBLE(thicknessS);
	if(_UC.thicknessB!=3.0) __OUT_DOUBLE(thicknessB);
	if(_UC.reconstruct!=1) __OUT_INT(reconstruct);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       VARIABLE COMPOSITION       *");
	__TITLE(line);
	g_free(line);
	if(_UC.firstGeneMax!=11) __OUT_INT(firstGeneMax);
	if(_UC._calctype_var){/*mandatory if varcomp*/
		__OUT_INT(minAt);
		__OUT_INT(maxAt);
	}
	if(_UC.fracTrans!=0.1) __OUT_DOUBLE(fracTrans);
	if(_UC.howManyTrans!=0.2) __OUT_DOUBLE(howManyTrans);
	if(_UC.specificTrans!=NULL) __OUT_BK_INT(specificTrans,"EndTransSpecific",_UC._nspetrans);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       METADYNAMICS       *");
	__TITLE(line);
	g_free(line);
	if(_UC.GaussianWidth!=0.) __OUT_DOUBLE(GaussianWidth);
	if(_UC.GaussianHeight!=0.) __OUT_DOUBLE(GaussianHeight);
	if(_UC.FullRelax!=2) __OUT_INT(FullRelax);
	if(_UC.maxVectorLength!=0.) __OUT_DOUBLE(maxVectorLength);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       PARTICLE SWARM OPTIMIZATION       *");
	__TITLE(line);
	g_free(line);
	if(_UC.PSO_softMut!=1.) __OUT_DOUBLE(PSO_softMut);
	if(_UC.PSO_BestStruc!=1.) __OUT_DOUBLE(PSO_BestStruc);
	if(_UC.PSO_BestEver!=1.) __OUT_DOUBLE(PSO_BestEver);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
is_w=0;
	line=g_strdup_printf("*       VCNEB       *");
	if(_UC.vcnebType!=110) __OUT_INT(vcnebType);
	if(_UC.numImages!=9) __OUT_INT(numImages);
	if(_UC.numSteps!=200) __OUT_INT(numSteps);
	if(_UC.optReadImages!=2) __OUT_INT(optReadImages);
	if((_UC._vcnebtype_method==1)&&(_UC.optimizerType!=1)) __OUT_INT(optimizerType);
	if((_UC._vcnebtype_method==2)&&(_UC.optimizerType!=2)) __OUT_INT(optimizerType);
	if(_UC.optRelaxType!=3) __OUT_INT(optRelaxType);
	if(_UC.dt!=0.05) __OUT_DOUBLE(dt);
	if(_UC.ConvThreshold!=0.003) __OUT_DOUBLE(ConvThreshold);
	if(_UC.VarPathLength!=0.) __OUT_DOUBLE(VarPathLength);
	if(_UC.K_min!=5.) __OUT_DOUBLE(K_min);
	if(_UC.K_max!=5.) __OUT_DOUBLE(K_max);
	if(_UC.Kconstant!=5.) __OUT_DOUBLE(Kconstant);
	if(_UC.optFreezing) __OUT_BOOL(optFreezing);
	if(_UC.optMethodCIDI!=0) __OUT_INT(optMethodCIDI);
	if(_UC.startCIDIStep!=100) __OUT_INT(startCIDIStep);
	if(_UC.pickupImages!=NULL) __OUT_BK_INT(pickupImages,"EndPickupImages",_UC._npickimg);
	if(_UC.FormatType!=2) __OUT_INT(FormatType);
	if(_UC.PrintStep!=1) __OUT_INT(PrintStep);
if(is_w==0) fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
else vfpos=ftell(vf);/* flag */
	line=g_strdup_printf("*       END OF FILE       *");
	__TITLE(line);
	g_free(line);
	fprintf(vf,"EOFEOFEOF\n");
	/*remove the garbage after EOFEOFEOF*/
	rewind(vf);
	line = file_read_line(vf);
	while(line){
		if(find_in_string("EOFEOFEOF",line) != NULL) break;
		fprintf(dest,"%s",line);
		g_free(line);
		line = file_read_line(vf);
	}
	if(line) g_free(line);
	/*there might be some garbage on the file after this point*/
	fclose(vf);
	fclose(dest);
	return 0;
#undef __TITLE
#undef __OUT_BOOL
#undef __OUT_INT
#undef __OUT_DOUBLE
#undef __OUT_STRING
#undef __OUT_BK_INT
#undef __OUT_BK_DOUBLE
#undef __OUT_BK_STRING
#undef __OUT_TMAT_DOUBLE
}
/*******************************/
/* Read Individual format file */
/*******************************/
gint read_individuals_uspex(gchar *filename, struct model_pak *model){
        FILE *vf;
	long int vfpos;
        gchar *line=NULL;
        uspex_output_struct *uspex_output=model->uspex;
        /* specific */
	gint idx=0;
	gchar *ptr;
	gint max_struct=0;
	gint red_index;
	/*since energy is sometimes eV, and sometimes eV/atom ?*/
	gboolean e_red=FALSE;
	/*start*/
        vf = fopen(filename, "rt");
        if (!vf) return 1;
	/*Have fitness? Better watch for each case...*/
	/*^^^^ I thought that calculationMethod = USPEX => fitness, but got defeated by EX19*/
	line = file_read_line(vf);
	if (find_in_string("Fitness",line) != NULL) _UO.have_fitness=TRUE;
	else _UO.have_fitness=FALSE;
	rewind(vf);
	/* deal with the problematic units */
	if(fetch_in_file(vf,"eV/atom")!=0) e_red=TRUE;
	rewind(vf);
        /* Skip header: FIX a _BUG_ when Individuals file sometimes has 2 line and sometimes only 1 */
        vfpos=ftell(vf);/* flag */
        line = file_read_line(vf);
        while(line){
                ptr=&(line[0]);
                while(*ptr==' ') ptr++;
                if(g_ascii_isdigit(*ptr)) {
                        fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
                        break;
                }
                vfpos=ftell(vf);/* flag */
                g_free(line);
                line = file_read_line(vf);
        }
        if(line==NULL) return -1;
	/* META calculation have a gen=0 Individual
	  unfortunately, it is *not* in the gather-
	 -POSCARS file, so we ignore gen=0 for now. */
	if(_UO.calc->calculationMethod==US_CM_META){
		g_free(line);
		line = file_read_line(vf);
		vfpos=ftell(vf);/* flag */
	}
	/*do a 1st pass to count stuctures*/
	_UO.num_struct=0;
        while (line){
		_UO.num_struct++;
		sscanf(line," %*i %i %*s",&idx);
		if(idx>max_struct) max_struct=idx;
		g_free(line);
		line = file_read_line(vf);
	}
	_UO.num_struct--;
	fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
	g_free(line);
	line = file_read_line(vf);
	/* now prepare arrays <- JOB=USPEX/300 example */
	if(_UO.num_struct==0) return 1;
	_UO.red_index=g_malloc(max_struct*sizeof(gint));
	for(idx=0;idx<max_struct;idx++) _UO.red_index[idx]=0;
	sscanf(line," %*[^[][%*[^]]] %lf %*f %*s",&(_UO.min_E));
	if(!e_red) _UO.min_E/=(gdouble)_UO.ind[0].natoms;
	_UO.max_E=_UO.min_E;
	idx=0;_UO.num_gen=1;
	while(line){
if(_UO.have_fitness){
	sscanf(line," %i %i%*[^[][%*[^]]] %lf %lf %*f %lf %*s",
		&(_UO.ind[idx].gen),&(red_index),&(_UO.ind[idx].energy),&(_UO.ind[idx].volume),&(_UO.ind[idx].fitness));
	if(isnan(_UO.ind[idx].fitness)) _UO.ind[idx].fitness=10000;
}else{
	sscanf(line," %i %i%*[^[][%*[^]]] %lf %lf %*s",
		&(_UO.ind[idx].gen),&(red_index),&(_UO.ind[idx].energy),&(_UO.ind[idx].volume));
}
		if(_UO.calc->calculationMethod==US_CM_META) red_index--;/*because of META "special" numbering*/
		if(_UO.ind[idx].gen<_UO.num_gen) _UO.ind[idx].gen++;
		if(_UO.ind[idx].gen>_UO.num_gen) _UO.num_gen=_UO.ind[idx].gen;
/*_VALGRIND_BUG_: _UO.red_index[red_index]=idx;*/
		_UO.red_index[red_index-1]=idx;
		if(!e_red) _UO.ind[idx].E=_UO.ind[idx].energy/_UO.ind[idx].natoms;
		else _UO.ind[idx].E=_UO.ind[idx].energy;
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX[Individual]: GEN=%i STRUCT=%i e=%lf E=%lf\n",_UO.ind[idx].gen,idx+1,_UO.ind[idx].energy,_UO.ind[idx].E);
#endif
		if(_UO.ind[idx].E<_UO.min_E) _UO.min_E=_UO.ind[idx].E;
		if(_UO.ind[idx].E>_UO.max_E) _UO.max_E=_UO.ind[idx].E;
		g_free(line);
		line = file_read_line(vf);
		idx++;
	}
	/*all done*/
	return 0;
}
/**********************************************************************************/
/* Read uspex output                                                              */
/* Require several files: OUTPUT.txt, Parameters.txt + calculation specific files.*/
/**********************************************************************************/
gint read_output_uspex(gchar *filename, struct model_pak *model){
	gchar *line;
	FILE *vf;
	long int vfpos;
	gchar *ptr;
	gchar *ptr2;
        gchar *res_folder;
	gchar *aux_file;
	gint probe;/*valgrind check*/
	gint ix;
	gint idx;
	gint jdx;
	gint min,max,med;
	gint job;
	/* results */
	FILE *f_src;
	FILE *f_dest;
	gdouble *e;
	gdouble *c;
	gdouble *tag;
	gdouble *tmp;
	gdouble compo;
	gint n_compo;
	gint skip=0;
	gint gen=1;
	gint num;
	gint natoms;
	gint species_index;
	gdouble min_E;
	gdouble min_F=0.;
	gdouble max_E;
	gchar *atoms;
	uspex_output_struct *uspex_output;
	uspex_calc_struct *uspex_calc;
	/* checks */
	g_return_val_if_fail(model != NULL, 1);
	g_return_val_if_fail(filename != NULL, 2);
	if (find_in_string("OUTPUT.txt",filename) == NULL){
		line = g_strdup_printf("ERROR: Please load USPEX OUTPUT.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		return -1;
	}
	/*allocs*/
	if(model->uspex!=NULL) g_free(model->uspex);
	model->uspex=g_malloc(sizeof(uspex_output_struct));
	uspex_output=model->uspex;
	/* TODO: check that selected file is OUTPUT.txt */
	vf = fopen(filename, "rt");
	if (!vf) return 1;
	error_table_clear();
/* --- get result folder*/
	res_folder=g_strdup (filename);
	ptr=g_strrstr (res_folder,"OUTPUT.txt");
	*ptr='\0';/*cut file at folder*/
#if DEBUG_USPEX_READ
	fprintf(stdout,"#DBG: USPEX_FOLDER=%s\n",res_folder);
#endif
/* --- setup environment */
	sysenv.render.show_energy = TRUE;
	
/* --- read the system OUTPUT.txt (will be confirmed by Parameters.txt)*/
	if(fetch_in_file(vf,"Version")==0){
		/*can't get version number (10.0.0 example will fail here)*/
		/* -> just produce a warning... file might actually be processed correctly*/
		line = g_strdup_printf("WARNING: USPEX version unsupported!\n");
		gui_text_show(WARNING, line);
		g_free(line);
	}else{
		rewind(vf);
		line = file_read_line(vf);
		while(line) {
			if (find_in_string("Version",line) != NULL){
				sscanf(line,"| Version %i.%i.%i %*s",&(max),&(med),&(min));
				_UO.version=min+10*med+100*max;
				if(_UO.version!=944){
					g_free(line);
					line = g_strdup_printf("WARNING: USPEX version %i unsupported!\n",_UO.version);
					gui_text_show(WARNING, line);

				}
				g_free(line);
				break;
			}
			g_free(line);
			line = file_read_line(vf);
		}
	}
	rewind(vf);
	job=0;
	if(fetch_in_file(vf,"Block for system description")==0) {
/*since VCNEB will fail here due to a non-unified OUTPUT.txt format
 *we need to be a little more permissive... */ 
		rewind(vf);
		if(fetch_in_file(vf,"VCNEB")==0) goto uspex_fail;/*definitely wrong*/
		job=-1;
	}
if(job>-1){
	/*next line contains Dimensionality*/
	ix=0;
	line = file_read_line(vf);
	/*unless from a previous version and it contains a lot of '-' character*/
	if (find_in_string("Dimensionality",line) == NULL) {
		g_free(line);
		line = file_read_line(vf);
	}
	sscanf(line," Dimensionality : %i ",&(ix));
	switch (ix){/*to be use later*/
		case -2:
		/*we do not deal with 2D crystal yet*/
		case 0:
		case 1:
		case 2:
		case 3:
		_UO.dim=ix;
		break;
		default:
		line = g_strdup_printf("ERROR: reading USPEX OUTPUT.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	job=ix*100;
	g_free(line);
	/*next line contains Molecular*/
	line = file_read_line(vf);
	sscanf(line," Molecular : %i %*s",&(ix));
	job+=(ix*10);
	if(ix==0) _UO.mol=FALSE;
	else if(ix==1) _UO.mol=TRUE;
	else{
		line = g_strdup_printf("ERROR: reading USPEX OUTPUT.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	g_free(line);
	/*next line contains Variable Composition*/
	line = file_read_line(vf);
	sscanf(line," Variable Composition : %i %*s",&(ix));
	job+=ix;
	if(ix==0) _UO.var=FALSE;
	else if(ix==1) _UO.var=TRUE;
	else{
		line = g_strdup_printf("ERROR: reading USPEX OUTPUT.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
#if DEBUG_USPEX_READ
	fprintf(stdout,"#DBG: USPEX OUTPUT calculationType=%i\n",job);
#endif
	g_free(line);
}/*in case of VCNEB, we don't have a correct job number*/
	line = file_read_line(vf);
	num=0;
	while(line){
		if(find_in_string("types of atoms in the system",line) != NULL) {
			/**/
			sscanf(line," There are %i %*s",&(num));
		}
		g_free(line);
		line = file_read_line(vf);
	}
	g_free(line);
/* --- and close */
	fclose(vf);vf=NULL;
/* --- READ Parameters.txt */
	aux_file = g_strdup_printf("%s%s",res_folder,"Parameters.txt");
	_UO.calc=read_uspex_parameters(aux_file,num);
	if(_UO.calc==NULL){
		line = g_strdup_printf("ERROR: reading USPEX REAL Parameter.txt file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		g_free(aux_file);
		goto uspex_fail;
	}
	g_free(aux_file);
	uspex_calc=_UO.calc;
	_UC.path=g_strdup(res_folder);
/* --- CHECK type from Parameters.txt matches that of OUTPUT.TXT */
	if((job!=-1)&&(job!=_UC.calculationType)){
/*since VCNEB will fail here due to a non-unified OUTPUT.txt format
 *we need to be a little more permissive... */
		line = g_strdup_printf("ERROR: inconsistent USPEX files (%i!=%i)!\n",job,_UC.calculationType);
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
if((_UC.calculationMethod==US_CM_USPEX)||(_UC.calculationMethod==US_CM_META)){
/* --- if calculationMethod={USPEX,META} */
	g_free(model->basename);/*FIX: _VALGRIND_BUG_*/
	model->basename=g_strdup_printf("uspex");
/* we have to open structure file first because of inconsistency in atom reporting in Individuals file*/
/* --- READ gatheredPOSCARS <- this is going to be the main model file */
/* ^^^ META is ill-defined: read first gatheredPOSCARS, then update with gatheredPOSCARS_relaxed in the future TODO*/
	aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
        vf = fopen(aux_file, "rt");
        if (!vf) {
                if(_UO.ind!=NULL) g_free(_UO.ind);
                line = g_strdup_printf("ERROR: can't open USPEX gatheredPOSCARS file!\n");
                gui_text_show(ERROR, line);
                g_free(line);
                goto uspex_fail;
        }
	/*count number of structures*/
	_UO.num_struct=0;
	line = file_read_line(vf);
	while(!feof(vf)){
		if((line[0]=='E')&&(line[1]=='A')) _UO.num_struct++;
		g_free(line);
		line = file_read_line(vf);
	}
	rewind(vf);
        model->num_frames=0;
        vfpos=ftell(vf);/* flag */
        line = file_read_line(vf);
	_UO.ind=g_malloc((_UO.num_struct)*sizeof(uspex_individual));
	idx=0;ix=-1;
        while(!feof(vf)){
                if((line[0]=='E')&&(line[1]=='A')) {
			idx=0;ix++;
                        fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
                        add_frame_offset(vf, model);
                        model->num_frames++;
                        g_free(line);
                        line = file_read_line(vf);
                }
		if(idx==5){/*calculate number of species and atoms*/
			_UO.ind[ix].atoms=g_malloc(_UC._nspecies*sizeof(gint));
			for(jdx=0;jdx<_UC._nspecies;jdx++) _UO.ind[ix].atoms[jdx]=0;/*init atoms*/
			_UO.ind[ix].natoms=0;
			/*get this structure nspecies*/
			ptr=&(line[0]);
			while(*ptr==' ') ptr++;/*skip initial space (if any)*/
			jdx=1;/*there is at least one species*/
			while(*ptr!='\0'){
				if(*ptr==' ') {
					jdx++;ptr++;
					while(*ptr==' ') ptr++;/*skip space*/
				}else ptr++;
			}
			if(jdx<_UC._nspecies){
				gchar *line2;/*double buffer*/
				/*we need to know which species is here*/
				line2 = file_read_line(vf);
				ptr=&(line[0]);
				ptr2=&(line2[0]);
				idx++;
				jdx=0;
				while(*ptr==' ') ptr++;/*skip initial space (if any)*/
				while(*ptr2==' ') ptr2++;/*skip initial space (if any)*/
				while((*ptr!='\0')&&(*ptr2!='\0')){
					/*get the correct jdx number*/
					jdx=0;/*we reset jdx to cope with 'OUT OF ORDER' poscar species, if any*/
					while((jdx<_UC._nspecies)&&(_UC.atomType[jdx]!=elem_symbol_test(ptr))) jdx++;
					_UO.ind[ix].atoms[jdx]=(gint)g_ascii_strtod(ptr2,NULL);
					_UO.ind[ix].natoms+=_UO.ind[ix].atoms[jdx];
					ptr++;ptr2++;
					while(g_ascii_isgraph(*ptr)) ptr++;/*go to next space/end*/
					while(g_ascii_isgraph(*ptr2)) ptr2++;/*go to next space/end*/
					ptr++;ptr2++;
					while(*ptr==' ') ptr++;/*skip space*/
					while(*ptr2==' ') ptr2++;/*skip space*/
				}
				g_free(line);
				line=line2;
			}else{
				/*get a new line (the number of each atoms)*/
				g_free(line);
				line = file_read_line(vf);
				idx++;
				ptr=&(line[0]);
				ptr2=ptr;jdx=0;
				do{/*get the number of each species (in the Parameters.txt order)*/
/*_VALGRIND_BUG_:			_UO.ind[ix].atoms[jdx]=g_ascii_strtod(ptr,&ptr2);*/
					probe=(gint)g_ascii_strtod(ptr,&ptr2);
					if(ptr2==ptr) break;
					_UO.ind[ix].atoms[jdx]=probe;
					_UO.ind[ix].natoms+=_UO.ind[ix].atoms[jdx];
					jdx++;
					ptr=ptr2;
				}while(1);
			}
		}
                vfpos=ftell(vf);/* flag */
		idx++;
                g_free(line);
                line = file_read_line(vf);
        }
	strcpy(model->filename,aux_file);// which means that we "forget" about OUTPUT.txt
        g_free(aux_file);
	fclose(vf);
/* --- READ Individuals <- information about all structures */
/* ^^^ META is ill-defined: read first Individuals, then update with Individual_relaxed and hope*/
	aux_file = g_strdup_printf("%s%s",res_folder,"Individuals");
	if(read_individuals_uspex(aux_file,model)!=0) {
		line = g_strdup_printf("ERROR: reading USPEX Individuals file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	g_free(aux_file);
	if(_UC.calculationMethod==US_CM_META) {
		aux_file = g_strdup_printf("%s%s",res_folder,"Individuals_relaxed");
		vf = fopen(aux_file, "rt");
		if (!vf) {
			g_free(aux_file);
		}else{
			/*we have an Individuals_relaxed, read it and update ONLY the structure defined inside*/
			vfpos=ftell(vf);/* flag */
			line = file_read_line(vf);
			while(line){
				ptr=&(line[0]);
				while(*ptr==' ') ptr++;
				if(g_ascii_isdigit(*ptr)) {
					fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
					break;
				}
				vfpos=ftell(vf);/* flag */
				g_free(line);
				line = file_read_line(vf);
			}
			g_free(line);
			fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
			line = file_read_line(vf);
			if(line==NULL) goto no_ind_relax;
			/*skip the gen=0 value!!*/
			g_free(line);line = file_read_line(vf);
			vfpos=ftell(vf);/* flag */
			/*no need to count, there will be less structure*/
			line = file_read_line(vf);
			while(line){
				sscanf(line," %*i %i %*s",&idx);
				if(_UO.have_fitness){
sscanf(line," %i %*i%*[^[][%*[^]]] %lf %lf %*f %lf %*s",
	&(_UO.ind[idx].gen),&(_UO.ind[idx].energy),&(_UO.ind[idx].volume),&(_UO.ind[idx].fitness));
if(isnan(_UO.ind[idx].fitness)) _UO.ind[idx].fitness=10000;
				}else{
sscanf(line," %i %*i%*[^[][%*[^]]] %lf %lf %*s",
	&(_UO.ind[idx].gen),&(_UO.ind[idx].energy),&(_UO.ind[idx].volume));
}
				g_free(line);
				line = file_read_line(vf);
			}
no_ind_relax:
			fclose(vf);
		}
	}
/* --- open the last frame */
/* ^^^ META is ill-defined: read first gatheredPOSCARS, then update with gatheredPOSCARS_relaxed in the future TODO*/
	aux_file = g_strdup_printf("%s%s",res_folder,"gatheredPOSCARS");
        vf = fopen(aux_file, "rt");g_free(aux_file);
        if (!vf) {/*very unlikely, we just opened it*/
                if(_UO.ind!=NULL) g_free(_UO.ind);
                line = g_strdup_printf("ERROR: can't open USPEX gatheredPOSCARS file!\n");
                gui_text_show(ERROR, line);
                g_free(line);
                goto uspex_fail;
        }
	read_frame_uspex(vf,model);/*open first frame*/
        fclose(vf);vf=NULL;
/* --- Prepare the summary graph */
/* ^^^ META is ill-defined: but we can read BESTIndividuals_relaxed instead of BESTIndividuals, hurrah!*/
	if(_UO.num_struct>_UO.num_gen){
		_UO.graph=graph_new("ALL", model);
		_UO.graph_best=graph_new("BEST", model);
		_UO.num_best=_UO.num_gen;/*I hope this one is always true*/
		_UO.best_ind=g_malloc((1+_UO.num_best)*sizeof(gint));
		for(ix=0;ix<(1+_UO.num_best);ix++) _UO.best_ind[ix]=0;
		/*NEW: read BESTIndividuals IDs first, when exists!*/
if(_UO.calc->calculationMethod==US_CM_META){
		aux_file = g_strdup_printf("%s%s",res_folder,"BESTIndividuals_relaxed");
		vf = fopen(aux_file, "rt");g_free(aux_file);
		if(!vf) {
			aux_file = g_strdup_printf("%s%s",res_folder,"BESTIndividuals");
			vf = fopen(aux_file, "rt");g_free(aux_file);
		}
}else{
		aux_file = g_strdup_printf("%s%s",res_folder,"BESTIndividuals");
		vf = fopen(aux_file, "rt");g_free(aux_file);
}
		probe=0;
if(vf){
		/*read BESTIndividuals*/
		/*skip header, avoid _BUG_ when Individual format have 1 OR 2 line header.*/
		vfpos=ftell(vf);/* flag */
		line = file_read_line(vf);
		while(line){
			ptr=&(line[0]);
			while(*ptr==' ') ptr++;
			if(g_ascii_isdigit(*ptr)) {
				fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
				break;
			}
			vfpos=ftell(vf);/* flag */
			g_free(line);
			line = file_read_line(vf);
		}
		if(line==NULL) goto no_best_ind;
	        /* META calculation have a gen=0 Individual
	          unfortunately, it is *not* in the gather-
	         -POSCARS file, so we ignore gen=0 for now. */
	        if(_UO.calc->calculationMethod==US_CM_META){
	                g_free(line);
	                line = file_read_line(vf);
	                vfpos=ftell(vf);/* flag */
	        }
		/*count best stuctures*/
		num=0;
		while (line){
			num++;
			g_free(line);
			line = file_read_line(vf);
		}
		num--;
		if(num>_UO.num_best) goto no_best_ind;/*num<_UO.num_best is allowed for unfinished calculation*/
		fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
		line = file_read_line(vf);
		/*now fill array*/
		while(line){
			sscanf(line," %i %i %*s",&gen,&ix);
			ix--;
			if(_UO.calc->calculationMethod==US_CM_META) ix--;/*META "special" numbering*/
			_UO.best_ind[gen]=ix;
#if DEBUG_USPEX_READ
        fprintf(stdout,"#DBG: USPEX GEN=%i BESTIndividuals=%i E=%f f=%f\n",
		gen,ix+1,_UO.ind[_UO.best_ind[gen]].energy,_UO.ind[_UO.best_ind[gen]].fitness);
#endif
			g_free(line);
			line = file_read_line(vf);
		}
		/*that's all <- we get values directly from Individuals*/
		probe=1;
no_best_ind:
		fclose(vf);
}
		skip=0;
		gen=1;
		ix=0;
		while(gen<=_UO.num_gen){
		  while(_UO.ind[ix].gen<gen) ix++;/*reach generation gen*/
		  num=0;
		  skip=ix;
if(probe==0){
		  min_E=_UO.ind[ix].E;
		  if(_UO.have_fitness) min_F=_UO.ind[ix].fitness;
		  _UO.best_ind[gen]=ix;
		  while(_UO.ind[ix].gen==gen){
			if(_UO.have_fitness){
				if(_UO.ind[ix].fitness<min_F){
					min_F=_UO.ind[ix].fitness;
					_UO.best_ind[gen]=ix;
				}
			}else{
				if(_UO.ind[ix].E<min_E){
					min_E=_UO.ind[ix].E;
					_UO.best_ind[gen]=ix;
				}
			}
			num++;
			ix++;
			if(ix>=_UO.num_struct) break;/*FIX _VALGRIND_BUG_*/
		  }
}else{
		  min_E=_UO.ind[_UO.best_ind[gen]].E;
		  if(_UO.have_fitness) min_F=_UO.ind[_UO.best_ind[gen]].fitness;
		  while(_UO.ind[ix].gen==gen){
			num++;
			ix++;
			if(ix>=_UO.num_struct) break;
		  }
}
		  e=g_malloc((1+num)*sizeof(gdouble));
		  ix=skip;
		  e[0]=(gdouble)num;
		  while(_UO.ind[ix].gen==gen){
			e[ix-skip+1]=_UO.ind[ix].E;
			ix++;
			if(ix>=_UO.num_struct) break;/*FIX _VALGRIND_BUG_*/
		  }
		  graph_add_borned_data(
			_UO.num_gen,e,0,_UO.num_gen,
			_UO.min_E-(_UO.max_E-_UO.min_E)*0.05,_UO.max_E+(_UO.max_E-_UO.min_E)*0.05,GRAPH_USPEX,_UO.graph);
		  gen++;
		  g_free(e);
		}
		/*prepare array*/
		e=g_malloc((1+_UO.num_gen)*sizeof(gdouble));
		e[0]=(gdouble)_UO.num_gen;
		max_E=_UO.ind[_UO.best_ind[1]].E;
		min_E=max_E;
		for(gen=0;gen<_UO.num_gen;gen++) {
			e[gen+1]=_UO.ind[_UO.best_ind[gen+1]].E;
			if(max_E<e[gen+1]) max_E=e[gen+1];
			if(min_E>e[gen+1]) min_E=e[gen+1];
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX_BEST: GEN=%i STRUCT=%i E=%lf e=%f\n",gen+1,1+_UO.best_ind[gen+1],e[gen+1],_UO.ind[_UO.best_ind[gen+1]].energy);
#endif
		}
		if(max_E==min_E) {
			min_E-=0.1;
			max_E+=0.1;
		}
		graph_add_borned_data(
			_UO.num_gen,e,0,_UO.num_gen,
			min_E-(max_E-min_E)*0.05,max_E+(max_E-min_E)*0.05,GRAPH_USPEX_BEST,_UO.graph_best);
		g_free(e);
	}
	/*set ticks*/
	if(_UO.num_gen>15) ix=5;
	else ix=_UO.num_gen+1;
	graph_set_xticks(TRUE,ix,_UO.graph);
	graph_set_yticks(TRUE,5,_UO.graph);
	graph_set_xticks(TRUE,ix,_UO.graph_best);
	graph_set_yticks(TRUE,5,_UO.graph_best);
/* --- variable composition */
if(_UO.var){
/* --- NEW: add a simple composition % graph for each species */
if(_UC._nspecies>1){
n_compo=2;
ix=0;
for(species_index=0;species_index<_UC._nspecies;species_index++){
	c=g_malloc(2*sizeof(gdouble));/*suppose at least 2 different compositions*/
	c[0]=(gdouble)_UO.ind[0].atoms[species_index] / (gdouble)_UO.ind[0].natoms;
	c[1]=c[0];
	compo=c[0];
	ix=0;
	while((compo==c[0])&&(ix<_UO.num_struct)){
		compo=(gdouble)_UO.ind[ix].atoms[species_index] / (gdouble)_UO.ind[ix].natoms;
		if(compo!=c[0]){
			if(compo>c[0]) c[1]=compo;
			else {
				c[1]=c[0];c[0]=compo;
			}
		}
		ix++;
	}
	if(c[1]==c[0]) continue;/*in VARCOMP calculation, one species can be kept fixed*/
	for(ix=0;ix<_UO.num_struct;ix++){
		compo=(gdouble)_UO.ind[ix].atoms[species_index] / (gdouble)_UO.ind[ix].natoms;
		if(compo<c[0]) {/*we have a new first element*/
			tmp=g_malloc((n_compo+1)*sizeof(gdouble));
				tmp=g_malloc((n_compo+1)*sizeof(gdouble));
				tmp[0]=compo;
				memcpy(&(tmp[1]),c,n_compo*sizeof(gdouble));
				g_free(c);
				c=tmp;
				n_compo++;
				continue;
			}
			if(compo>c[n_compo-1]){/*we have a new last element*/
				tmp=g_malloc((n_compo+1)*sizeof(gdouble));
				tmp[n_compo]=compo;
				memcpy(&(tmp[0]),c,n_compo*sizeof(gdouble));
				g_free(c);
				c=tmp;
				n_compo++;
				continue;
			}
			for(jdx=1;jdx<n_compo;jdx++){
				if((compo>c[jdx-1])&&(compo<c[jdx])){
				/*insert a new element at position jdx*/
					tmp=g_malloc((n_compo+1)*sizeof(gdouble));
					tmp[jdx]=compo;
					memcpy(&(tmp[0]),c,(jdx)*sizeof(gdouble));
					memcpy(&(tmp[jdx+1]),&(c[jdx]),(n_compo-jdx)*sizeof(gdouble));
					g_free(c);
					c=tmp;
					n_compo++;
					jdx=n_compo;/*exit loop*/
					continue;
				}
				
			}
		}
		/*build the energy array*/
		e=g_malloc(n_compo*sizeof(gdouble));
		tag=g_malloc(n_compo*sizeof(gdouble));
		max_E=_UO.min_E;
		for(ix=0;ix<n_compo;ix++) e[ix]=_UO.max_E;
		for(ix=0;ix<_UO.num_struct;ix++){
			for(jdx=0;jdx<n_compo;jdx++){
				compo=(gdouble)_UO.ind[ix].atoms[species_index] / (gdouble)_UO.ind[ix].natoms;
				if(c[jdx]==compo){
					if(_UO.ind[ix].E < e[jdx]) {
						e[jdx] = _UO.ind[ix].E;
						tag[jdx] = (gdouble) ix;
						if(e[jdx]>max_E) max_E=e[jdx];
					}
				}
			}
		}
		/* prepare composition diagram*/
		line=g_strdup_printf("COMP_%s",elements[_UC.atomType[species_index]].symbol);
		_UO.graph_comp=graph_new(line, model);
                graph_add_borned_data(
			n_compo,tag,0,1,_UO.min_E-(max_E-_UO.min_E)*0.05,max_E+(max_E-_UO.min_E)*0.05,GRAPH_USPEX_2D,_UO.graph_comp);
		graph_add_borned_data(
			n_compo,c,0,1,_UO.min_E-(max_E-_UO.min_E)*0.05,max_E+(max_E-_UO.min_E)*0.05,GRAPH_USPEX_2D,_UO.graph_comp);
		graph_add_borned_data(
			n_compo,e,0,1,_UO.min_E-(max_E-_UO.min_E)*0.05,max_E+(max_E-_UO.min_E)*0.05,GRAPH_USPEX_2D,_UO.graph_comp);
                g_free(line);
		g_free(tag);
		g_free(c);
		g_free(e);
        	/*set ticks*/
        	graph_set_xticks(TRUE,3,_UO.graph_comp);
        	graph_set_yticks(TRUE,5,_UO.graph_comp);

	}

}/*variable composition with only 1 species...*/
}/*not a VARCOMP calculation*/

	/*refresh model*/
	tree_model_refresh(model);
        model->redraw = TRUE;
        /* always show this information */
        gui_text_show(ITALIC,g_strdup_printf("USPEX: %i structures detected.\n",model->num_frames));
        model_prep(model);


}else if(_UC.calculationMethod==US_CM_VCNEB){
/* --- tentative at adding VCNEB, which format is different**/
/* ^^^ *BUT unfortunately, we NEED a vasp5 output (VCNEB stayed at VASP4) */
/*      thus we create a new file called transitionPath_POSCARs5 which is a
 *      translation of transitionPath_POSCARs for VASP5...*/
	atoms = g_strdup_printf("%s",elements[_UC.atomType[0]].symbol);
	for(idx=1;idx<_UC._nspecies;idx++){
		line = g_strdup_printf("%s %s",atoms,elements[_UC.atomType[idx]].symbol);
		g_free(atoms);
		atoms=line;
	}
/* open transitionPath_POSCARs and start to convert to transitionPath_POSCARs5 */
	aux_file = g_strdup_printf("%s%s",res_folder,"transitionPath_POSCARs");
	f_src = fopen(aux_file, "rt");
	if(!f_src) {
		line = g_strdup_printf("ERROR: can't open USPEX transitionPath_POSCARs file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	g_free(aux_file);
	aux_file = g_strdup_printf("%s%s",res_folder,"transitionPath_POSCARs5");
	f_dest = fopen(aux_file, "w");
	if(!f_dest) {
		line = g_strdup_printf("ERROR: can't WRITE into USPEX transitionPath_POSCARs5 file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
/*get the total number of atoms*/
	for(idx=0;idx<5;idx++){
		line = file_read_line(f_src);
		g_free(line);
	}
	natoms=1;/*there is at least one atom*/
	ptr=&(line[0]);
	do{
		natoms+=(gint)g_ascii_strtod(ptr,&ptr2);
		if(ptr2==ptr) break;
		ptr=ptr2;
	}while(1);
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX[VCNEB]: NATOMS=%i ATOMS=%s\n",natoms,atoms);
#endif
	/*translate input VASP4 into proper VASP5*/
	rewind(f_src);
	idx=0;
	_UO.num_struct=0;
	line = file_read_line(f_src);
	while(line){
		fprintf(f_dest,"%s",line);
		if (find_in_string("Image",line) != NULL) {
			idx=0;
			_UO.num_struct++;
		}
		if(idx==4) fprintf(f_dest,"%s\n",atoms);
		idx++;
		g_free(line);
		line = file_read_line(f_src);
	}
	fclose(f_src);
	fclose(f_dest);
	g_free(atoms);
	_UO.ind=g_malloc((_UO.num_struct)*sizeof(uspex_individual));
/* --- read energies from the Energy file */
	g_free(aux_file);
	aux_file = g_strdup_printf("%s%s",res_folder,"Energy");
	vf=fopen(aux_file,"rt");
	if (!vf) {
		g_free(_UO.ind);
		line = g_strdup_printf("ERROR: can't open USPEX Energy file!\n");
		gui_text_show(ERROR, line);
		g_free(line);
		goto uspex_fail;
	}
	/*first get to the first data line*/
	line = file_read_line(vf);
	ptr=&(line[0]);
	while((ptr)&&(*ptr==' ')) ptr++;/*skip blanks*/
	while(!g_ascii_isdigit(*ptr)){
		g_free(line);
		line = file_read_line(vf);
		ptr=&(line[0]);
		while((ptr)&&(*ptr==' ')) ptr++;/*skip blanks*/
	}/*line should now point to the first data line*/
	/*get first energy*/
	_UO.best_ind=g_malloc((1+_UO.num_struct)*sizeof(gint));
	e = g_malloc((_UO.num_struct+1)*sizeof(gdouble));
	e[0]=(gdouble) _UO.num_struct;/*first element is the size*/
	sscanf(line," %*i %lf %*s",&(e[1]));
	e[1] /= (gdouble)natoms;
	_UO.min_E=e[1];
	_UO.max_E=e[1];
	for(idx=0;idx<_UO.num_struct;idx++){
		if(!line) break;/*<- useful? */
		sscanf(line," %*i %lf %*s",&(_UO.ind[idx].energy));
		_UO.ind[idx].atoms=g_malloc(_UC._nspecies*sizeof(gint));
		_UO.ind[idx].natoms=natoms;/*<-constant*/
		_UO.ind[idx].E = _UO.ind[idx].energy / (gdouble)_UO.ind[idx].natoms;
		_UO.ind[idx].gen=idx;/*<- this is actually *not* true*/
		_UO.best_ind[idx+1]=idx;
		e[idx+1]=_UO.ind[idx].E;
		if(e[idx+1]<_UO.min_E) _UO.min_E=e[idx+1];
		if(e[idx+1]>_UO.max_E) _UO.max_E=e[idx+1];
#if DEBUG_USPEX_READ
fprintf(stdout,"#DBG: USPEX[VCNEB] Image %i: natoms=%i energy=%lf e=%lf\n",_UO.ind[idx].gen,_UO.ind[idx].natoms,_UO.ind[idx].energy,e[idx]);
#endif
		g_free(line);
		line = file_read_line(vf);
	}
	fclose(vf);vf=NULL;
	g_free(aux_file);
	/*create the reduced index table*/
	_UO.red_index=g_malloc(_UO.num_struct*sizeof(gint));
	for(idx=0;idx<_UO.num_struct;idx++) _UO.red_index[idx]=idx;
	aux_file = g_strdup_printf("%s%s",res_folder,"transitionPath_POSCARs5");
/* --- reopen the newly created transitionPath_POSCARs5 file for reading <- this will be our new model file*/
	vf=fopen(aux_file,"rt");
        if (!vf) {/*very unlikely: we just create it*/
		g_free(_UO.ind);
		line = g_strdup_printf("ERROR: can't open USPEX transitionPath_POSCARs5 file!\n");
                gui_text_show(ERROR, line);
                g_free(line);
                goto uspex_fail;
        }
        model->num_frames=0;
        vfpos=ftell(vf);/* flag */
        line = file_read_line(vf);
	idx=0;
        while(!feof(vf)){
		if (find_in_string("Image",line) != NULL) {
			idx=0;
                        fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
                        add_frame_offset(vf, model);
                        model->num_frames++;
                        g_free(line);
                        line = file_read_line(vf);
		}
                vfpos=ftell(vf);/* flag */
                g_free(line);
                line = file_read_line(vf);
		idx++;
        }
        strcpy(model->filename,aux_file);// which means that we "forget" about OUTPUT.txt
        g_free(aux_file);
        rewind(vf);
        read_frame_uspex(vf,model);/*open first frame*/
        fclose(vf);vf=NULL;

	/*do a "best" graph with each image*/
        _UO.graph_best=graph_new("PATH", model);
	_UO.num_gen=_UO.num_struct;
        graph_add_borned_data(
                _UO.num_gen,e,0,_UO.num_gen,
                _UO.min_E-(_UO.max_E-_UO.min_E)*0.05,_UO.max_E+(_UO.max_E-_UO.min_E)*0.05,GRAPH_USPEX_BEST,_UO.graph_best);
        g_free(e);
        /*set ticks*/
        if(_UO.num_gen>15) ix=5;
        else ix=_UO.num_gen+1;
        graph_set_xticks(TRUE,ix,_UO.graph_best);
        graph_set_yticks(TRUE,5,_UO.graph_best);


        /*refresh model*/
        tree_model_refresh(model);
        model->redraw = TRUE;
        /* always show this information */
        gui_text_show(ITALIC,g_strdup_printf("USPEX: %i structures detected.\n",model->num_frames));
        model_prep(model);



}else{/*calculation method is not supported*/
	line = g_strdup_printf("ERROR: USPEX support is limited to calculationMethod={USPEX,META,VCNEB} for now!\n");
	gui_text_show(ERROR, line);
	g_free(line);
	goto uspex_fail;
}
	/*g_free some data*/
	g_free(res_folder);
	/*update some properties*/
	model->num_species=_UO.calc->_nspecies;

	/*end*/
	error_table_print_all();
	return 0;
uspex_fail:
	line = g_strdup_printf("ERROR: loading USPEX failed!\n");
	gui_text_show(ERROR, line);
	g_free(line);
	return 3;
}
gint read_frame_uspex(FILE *vf, struct model_pak *model){
	gchar *line;
	uspex_output_struct *uspex_output=model->uspex;
	long int vfpos;/*to counter _BUG_ where loading a raw frame does not update model->cur_frame*/
	gchar *ptr;
	gint idx=0;
	/*read frame*/
        g_assert(vf != NULL);

vfpos=ftell(vf);/* flag */
/*_BUG_ model->cur_frame not updated, use frame title instead*/
line = file_read_line(vf);
if(line==NULL) return -1;
ptr=&(line[0]);
while(g_ascii_isalpha(*ptr)) ptr++;
sscanf(ptr,"%i%*s",&idx);
g_free(line);/*FIX: _VALGRIND_BUG_*/
fseek(vf,vfpos,SEEK_SET);/* rewind to flag */
        if(vasp_load_poscar5(vf,model)<0) return 3;
/* in case of meta idx is *not* the real number of the structure */
/* now the index should always be _UO.red_index[idx-1] <- I THINK*/
	line=g_strdup_printf("%lf eV",_UO.ind[_UO.red_index[idx-1]].energy);
	property_add_ranked(3, "Energy", line, model);
	model->volume=_UO.ind[_UO.red_index[idx-1]].volume;
	g_free(line);line=g_strdup_printf("%lf uV",model->volume);
	property_add_ranked(4, "Volume", line, model);
	g_free(line);
	/*note: Formula property rank is 7*/
	line=g_strdup_printf("%i",idx);
	property_add_ranked(8, "Structure", line, model);
	g_free(line);
if(_UO.have_fitness){
	line=g_strdup_printf("%f f",_UO.ind[_UO.red_index[idx-1]].fitness);
	property_add_ranked(9, "Fitness", line, model);
	g_free(line);
}
	return 0;
}

#undef _UO

