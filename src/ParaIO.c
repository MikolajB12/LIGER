/*  Liger - Calculates the distribution of galaxies (or any targets) on the backward light cone of any observer taking into account leading-order GR effects in the weak-field regime
    Copyright (C) 2017  Mikolaj Borzyszkowski, Argelander Institut für Astronomie, University of Bonn

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

//  This file contains the routine which reads the parameters from the Parameter file.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./vardef.h"

#define NUMPARAS 20
#define TEXTPARA 0
#define FLOATPARA 1
#define INTPARA 2
#define TRUTHPARA 3
#define FLOATARRAYPARA 4

int read_ParaFile(char *fname){
char Paras[NUMPARAS][100],ReadLine[500],ReadParaBuf[200],ReadVal1[200],ReadVal2[200],ReadVal3[200];
FILE *Pfile;
int ParaType[NUMPARAS],Npara,i,Nread;
void *ParaAdr[NUMPARAS];

Npara=0;

sprintf(Paras[Npara],"DirOutput");
ParaType[Npara]=TEXTPARA;
ParaAdr[Npara]=run.outpath;
++Npara;

sprintf(Paras[Npara],"FileNameOutput");
ParaType[Npara]=TEXTPARA;
ParaAdr[Npara]=run.outname;
++Npara;

sprintf(Paras[Npara],"DirGravPotential");
ParaType[Npara]=TEXTPARA;
ParaAdr[Npara]=run.gravpath;
++Npara;

sprintf(Paras[Npara],"DirSimulation");
ParaType[Npara]=TEXTPARA;
ParaAdr[Npara]=run.simpath;
++Npara;

sprintf(Paras[Npara],"NumSimSnapshots");
ParaType[Npara]=INTPARA;
ParaAdr[Npara]=&run.NSnap;
++Npara;

sprintf(Paras[Npara],"NumThreads");
ParaType[Npara]=INTPARA;
ParaAdr[Npara]=&run.Nthreads;
++Npara;

sprintf(Paras[Npara],"NGravGrid");
ParaType[Npara]=INTPARA;
ParaAdr[Npara]=&run.Ngrid;
++Npara;

sprintf(Paras[Npara],"IfWritePotential");
ParaType[Npara]=TRUTHPARA;
ParaAdr[Npara]=&run.ifOutputPotential;
++Npara;

sprintf(Paras[Npara],"NTargetsInMem");
ParaType[Npara]=INTPARA;
ParaAdr[Npara]=&run.NPartInMem;
++Npara;

sprintf(Paras[Npara],"MemBufFactor");
ParaType[Npara]=FLOATPARA;
ParaAdr[Npara]=&run.OutBuf;
++Npara;

sprintf(Paras[Npara],"Maxdist");
ParaType[Npara]=FLOATPARA;
ParaAdr[Npara]=&run.maxdist;
++Npara;

sprintf(Paras[Npara],"Mindist");
ParaType[Npara]=FLOATPARA;
ParaAdr[Npara]=&run.mindist;
++Npara;

sprintf(Paras[Npara],"ObsPos");
ParaType[Npara]=FLOATARRAYPARA;
ParaAdr[Npara]=&run.ObsPos[0];
++Npara;

if((Pfile=fopen(fname, "r"))) {      
	while (fgets(ReadLine, sizeof(ReadLine),Pfile) != NULL){
		Nread=sscanf(ReadLine,"%s%s%s%s",ReadParaBuf,ReadVal1,ReadVal2,ReadVal3);
		if (Nread < 2) continue;
		if (ReadParaBuf[0]==*"%") continue;

		for (i=0;i<Npara;++i) if (strcmp(ReadParaBuf, Paras[i]) == 0) break;
		if (i==Npara) {
			printf("	Ignoring unkown parameter %s\n",ReadParaBuf);
			continue;
		}
		switch (ParaType[i]) {
			case FLOATPARA:
			 ParaType[i]=-1;
			 if (sscanf(ReadVal1,"%lf",((double *)ParaAdr[i]))!=1){
				printf("	Error reading %s, should be an floating-point value!\n",ReadParaBuf);
				exit(0);
			 } break;

			case FLOATARRAYPARA:
			 if (Nread < 4) {
				printf("	Expect more numerical values for %s\n",ReadParaBuf);
				exit(0);
			 }
			 ParaType[i]=-1;
			 if (sscanf(ReadVal1,"%lf",&(((double *)ParaAdr[i])[0]))!=1){
				printf("	Error reading %s, should be an floating-point value!\n",ReadParaBuf);
				exit(0);
			 }
			 if (sscanf(ReadVal2,"%lf",&(((double *)ParaAdr[i])[1]))!=1){
				printf("	Error reading %s, should be an floating-point value!\n",ReadParaBuf);
				exit(0);
			 }
			 if (sscanf(ReadVal3,"%lf",&(((double *)ParaAdr[i])[2]))!=1){
				printf("	Error reading %s, should be an floating-point value!\n",ReadParaBuf);
				exit(0);
			 } break;

			case TEXTPARA:
			 ParaType[i]=-1;
			 strcpy(ParaAdr[i], ReadVal1); break;

			case TRUTHPARA:
			 ParaType[i]=-1;
			 if (strcmp(ReadVal1, "True") == 0) ((int *)ParaAdr[i])[0]=1;
			 else if (strcmp(ReadVal1, "False") == 0) ((int *)ParaAdr[i])[0]=0;
			 else {
				printf("	Only a accept True or False for %s\n",ReadParaBuf);
				exit(0);
			 } break;

			case INTPARA:
			 ParaType[i]=-1;
			 if (sscanf(ReadVal1,"%i",((int *)ParaAdr[i]))!=1){
				printf("	Error reading %s, should be an integer value!\n",ReadParaBuf);
				exit(0);
			 } break;

			case -1:
			 printf("	Parameter %s appears more then once!\n",ReadParaBuf);
			 exit(0); break;
		}
	}
	fclose(Pfile);
} else {
	printf("	Cannot read file with parameters from: %s\n",fname);
	exit(0);
}

Nread=0;
for (i=0;i<Npara;++i) {
	if (ParaType[i]!=-1) printf("	Missing parameter: %s\n",Paras[i]);
	else ++Nread;
}
if (Nread!=Npara) {
	printf("%i parameters missing.\n",Npara-Nread);
	exit(0);
}

return Nread;
}

#undef NUMPARAS
#undef TEXTPARA
#undef FLOATPARA
#undef INTPARA
#undef TRUTHPARA
#undef FLOATARRAYPARA
