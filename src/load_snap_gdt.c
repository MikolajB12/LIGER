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

//  This file contains the routines which read the simulation data and the list of galaxies.
//  These are adapted versions of the procedures published along with GADGET2, Springel, V. (2005)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./vardef.h"
#include "./load_snap_gdt.h"

struct header_snapshot load_header_gdt(char *fname){
	struct header_snapshot head;
	FILE *fd;
	char blkname[4];
	int dummy,SnapFormat=1;

	if(!(fd = fopen(fname, "r"))){
		printf("can't open file `%s`\n", fname);
		exit(0);
	}

	printf("	reading header from `%s' ...", fname);
	fflush(stdout);

	fread(&dummy, sizeof(dummy), 1, fd);
	if (dummy==8) SnapFormat=2; else if (dummy==256) SnapFormat=1; else {printf("Unknown Snapshot format!\nFirst check sum=%i\n",dummy); fclose(fd); exit(0);}
	if (SnapFormat==2){
		fread(blkname,sizeof(char),4,fd);
		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);
	}
	fread(&head, sizeof(head), 1, fd);
	fread(&dummy, sizeof(dummy), 1, fd);
	fclose(fd);

	//printf("Header:\nNPart:		%i	%i	%i	%i	%i	%i\n",header1.npart[0],header1.npart[1],header1.npart[2],header1.npart[3],header1.npart[4],header1.npart[5]);
	//printf("NPartTot:	%i	%i	%i	%i	%i	%i\n",header1.npartTotal[0],header1.npartTotal[1],header1.npartTotal[2],header1.npartTotal[3],header1.npartTotal[4],header1.npartTotal[5]);
	//printf("Mass:		%e	%e	%e	%e	%e	%e\n",header1.mass[0],header1.mass[1],header1.mass[2],header1.mass[3],header1.mass[4],header1.mass[5]);
	//printf("Time: %f	Red: %f	SFR: %i	FB: %i	Cool: %i	Nfiles: %i\n",header1.time,header1.redshift,header1.flag_sfr,header1.flag_feedback,header1.flag_cooling,header1.num_files);
	//printf("Box: %e	OM: %e	OL: %e	H: %e\n",header1.BoxSize,header1.Omega0,header1.OmegaLambda,header1.HubbleParam);

	printf("done.\n");
	return head;
}



// Here a gadget simulation snapshot file is read into memory.
int load_snapshot_gdt(char *fname, int files,struct particle_data **Pptr)
{
struct header_snapshot header1;
FILE *fd;
char buf[200],blkname[4];
int i, k, dummy, ntot_withmasses,SnapFormat;
int n, pc, pc_new, NumPart, Ngas;
//int j, t, n, off, pc_sph;
struct particle_data *P;

SnapFormat=1;

for(i = 0, pc = 0; i < files; i++, pc = pc_new) {
	if(files > 1)
		sprintf(buf, "%s.%d", fname, i);
	else
		sprintf(buf, "%s", fname);

	if(!(fd = fopen(buf, "r"))){
		printf("can't open file `%s`\n", buf);
		exit(0);
	}

	printf("	reading `%s' ...", buf);
	fflush(stdout);

	fread(&dummy, sizeof(dummy), 1, fd);
	if (dummy==8) SnapFormat=2; else if (dummy==256) SnapFormat=1; else {printf("Unknown Snapshot format!\nFirst check sum=%i\n",dummy); fclose(fd); exit(0);}
	if (SnapFormat==2){
		fread(blkname,sizeof(char),4,fd);
		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);
	}
	fread(&header1, sizeof(header1), 1, fd);
	fread(&dummy, sizeof(dummy), 1, fd);
	//printf("%i\n",dummy);

	//printf("Header:\nNPart:		%i	%i	%i	%i	%i	%i\n",header1.npart[0],header1.npart[1],header1.npart[2],header1.npart[3],header1.npart[4],header1.npart[5]);
	//printf("NPartTot:	%i	%i	%i	%i	%i	%i\n",header1.npartTotal[0],header1.npartTotal[1],header1.npartTotal[2],header1.npartTotal[3],header1.npartTotal[4],header1.npartTotal[5]);
	//printf("Mass:		%e	%e	%e	%e	%e	%e\n",header1.mass[0],header1.mass[1],header1.mass[2],header1.mass[3],header1.mass[4],header1.mass[5]);
	//printf("Time: %f	Red: %f	SFR: %i	FB: %i	Cool: %i	Nfiles: %i\n",header1.time,header1.redshift,header1.flag_sfr,header1.flag_feedback,header1.flag_cooling,header1.num_files);
	//printf("Box: %e	OM: %e	OL: %e	H: %e\n",header1.BoxSize,header1.Omega0,header1.OmegaLambda,header1.HubbleParam);

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

	if(files == 1) {
		for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
			NumPart += header1.npart[k];
		Ngas = header1.npart[0];
	} else {
		for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
			NumPart += header1.npartTotal[k];
		Ngas = header1.npartTotal[0];
	}

	for(k = 0, ntot_withmasses = 0; k < 6; k++) {
		if(header1.mass[k] == 0) ntot_withmasses += header1.npart[k];
	}

	if(i == 0) {
		allocate_memory(NumPart,Pptr);
		P=Pptr[0];
	}

	SKIP;
	if (SnapFormat==2) {
		fread(&blkname[0],sizeof(char),4,fd);
		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);
	}
	for(k = 0, pc_new = pc; k < 6; k++) {
		for(n = 0; n < header1.npart[k]; n++) {
			fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
			pc_new++;
		}
	}
	SKIP;
//Following we could read the velocity, particle ID's or 
      fclose(fd);
    }

  printf("done.\n");
  fflush(stdout);
}


// Here we read a gadget simulation snapshot partially (currently only files=1 is working).
int load_snapshot_gdt_partial(char *fname,int files,int NStart,int NSelect,struct particle_data **Pptr)
{
FILE *fd;
char buf[200];
double SqrtExpan;
int dummy;
int i, j, k, ntot_withmasses,Nfile;
int t, n, off, pc, pc_new, pc_sph, NumPart, Ngas;
struct header_snapshot header1;
struct particle_data *P=Pptr[0];

for(i = 0, pc = 0; i < files; i++, pc = pc_new) {
	if(files > 1) sprintf(buf, "%s.%d", fname, i);
	else sprintf(buf, "%s", fname);

	if(!(fd = fopen(buf, "r"))) {
		return -1;
		//printf("can't open file `%s`\n", buf);
		//exit(0);
	}

	fread(&dummy, sizeof(dummy), 1, fd);
	fread(&header1, sizeof(header1), 1, fd);
	fread(&dummy, sizeof(dummy), 1, fd);

	#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

	if(files == 1) {
		for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
			NumPart += header1.npart[k];
		Ngas = header1.npart[0];
	} else {
		for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
			NumPart += header1.npartTotal[k];
		Ngas = header1.npartTotal[0];
	}

	for(k = 0, Nfile = 0; k < 6; k++) Nfile += header1.npart[k];

	for(k = 0, ntot_withmasses = 0; k < 6; k++) {
		if(header1.mass[k] == 0) ntot_withmasses += header1.npart[k];
	}

	if(i == 0) {
		allocate_memory_partial(NSelect,Pptr);
		P=Pptr[0];
	}

	if ((pc+Nfile<NStart)||(pc>NStart+NSelect)) {
		for (k=0;k<6;k++) pc+=header1.npart[k];
		fclose(fd);
		continue;
	}

	SKIP;
	pc_new = pc;
	n=0;
	if (NStart>pc) {
		fseek(fd,3*sizeof(float)*(NStart-pc),SEEK_CUR);
		pc_new +=NStart-pc;
		n=NStart-pc;
	}
	for(; (n < Nfile)&&(pc_new < NStart+NSelect); n++) { //n runs over the number of particles in the file
		fread(&P[pc_new-NStart].Pos[0], sizeof(float), 3, fd); //pc_new is total particle count
		pc_new++;
	}
	if (n < Nfile) {
		fseek(fd,3*sizeof(float)*(Nfile-n),SEEK_CUR);
		pc_new +=Nfile-n;
	}
	SKIP;

	SKIP;
	pc_new = pc;
	n=0;
	if (NStart>pc) {
		fseek(fd,3*sizeof(float)*(NStart-pc),SEEK_CUR);
		pc_new +=NStart-pc;
		n=NStart-pc;
	}
	for(; (n < Nfile)&&(pc_new < NStart+NSelect); n++) { //n runs over the number of particles in the file
		fread(&P[pc_new-NStart].Vel[0], sizeof(float), 3, fd); //pc_new is total particle count
		pc_new++;
	}
	if (n < Nfile) {
		fseek(fd,3*sizeof(float)*(Nfile-n),SEEK_CUR);
		pc_new +=Nfile-n;
	}
	SKIP;

      fclose(fd);
    }

  SqrtExpan = sqrt(header1.time);
  for (n=0;n<NSelect;n++) for (k=0;k<3;k++) P[n].Vel[k]*=SqrtExpan;

  return 0;
}




// this routine allocates the memory for the particle data.
int allocate_memory(int NumPart,struct particle_data **P)
{

  if(!(P[0] = malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

}

int allocate_memory_partial(int NSelection,struct particle_data **P)
{

  if(!(P[0] = malloc(NSelection * sizeof(struct particle_data))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }
}


