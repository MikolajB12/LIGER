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
//  These routines only redirect their input, this is doe to facilitate the use of custom routines to handle differet format of input data.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "./vardef.h"
#include "./load_snap_gdt.h"

struct header_snapshot load_snapshot_header(int FileNr){
char fname[MAXPATHLENGTH];
sprintf(fname,"%s%03d",run.simpath,FileNr);
return load_header_gdt(fname);
}

struct header_snapshot load_galaxies_header(int FileNr){
char fname[MAXPATHLENGTH];
sprintf(fname,"%s%03d",run.simpath,FileNr);
return load_header_gdt(fname);
}

int load_snapshot(int FileNr, int files,struct particle_data **Pptr){
char fname[MAXPATHLENGTH];
sprintf(fname,"%s%03d",run.simpath,FileNr);
#ifdef MEASURETIME
 IOTime[0]-=time(NULL);
 int result=load_snapshot_gdt(fname,files,Pptr);
 IOTime[0]+=time(NULL);
 return result;
#else
 return load_snapshot_gdt(fname,files,Pptr);
#endif
}

//Currently for Gadget only files=1 is possible.
int load_galaxies_partial(int FileNr,int files,int NStart,int NSelect,struct particle_data **Pptr){
char fname[MAXPATHLENGTH];
sprintf(fname,"%s%03d",run.simpath,FileNr);
#ifdef MEASURETIME
 #pragma omp atomic
  IOTime[1]-=time(NULL);
 int result=load_snapshot_gdt_partial(fname,files,NStart,NSelect,Pptr);
 #pragma omp atomic
  IOTime[1]+=time(NULL);
 return result;
#else
 return load_snapshot_gdt_partial(fname,files,NStart,NSelect,Pptr);
#endif
}

