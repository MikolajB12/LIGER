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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef NOOPENMP
 #include <omp.h>
#endif
#include <fftw3.h> 
#include <unistd.h>
#include "./src/vardef.h"
#include "./src/ParaIO.h"
#include "./src/load_snap.h"
#include "./src/integration.h"
#include "./src/potential.h"
#include "./src/light_rays.h"

int main(int argc, char **argv)
{
	char path[MAXPATHLENGTH],gravpath[MAXPATHLENGTH];
	int a,b,c,d,i,j,bStart,bStop,Ngrid,Ngrid3,NMaxOutput,NextPart=0,RealStepPart,NumGal;
	double gridsize,OmM,OmL,boxsize;
	double *rangemin,*rangemax,*Hubarray;
	float *PotHold,LCsize;
	struct potential_container **Pot;
	struct particle_data **AllP;
	struct header_snapshot *Head,Header;
	struct light_rays_path *ray;
	struct light_cone_object *AllConePos,ConePos;
	FILE *out;

	#ifdef MEASURETIME
	 TotalTime=-time(NULL);
	 for (i=0;i<2;++i) PartTime[i]=0;
	 for (i=0;i<2;++i) IOTime[i]=0;
	#endif

	printf("\n");
	printf("       L                 ______     \n");
	printf("             G  _       /      \\    \n");
	printf(" <__________   / \\ R   /  v1.0  \\   \n");
	printf("            \\_/ E \\   /          \\  \n");
	printf("          I        \\_/            \\ \n");
	printf("\n");


	if (argc < 2) {
		printf("Please supply a paramter file!\nOtherwise I don't know what to do.\n");
		exit(0);
	}
	printf("\n.....................................................\n");
	printf("Reading %s\n",argv[1]);
	printf("	Found %i parameters.\n",read_ParaFile(argv[1]));
	fflush(stdout);

	#ifndef NOOPENMP
	 if (0==fftw_init_threads()) {
		printf("Failed to initialize fftw threads! Let's try without.\n");
		//exit(0);
	 } else fftw_plan_with_nthreads(run.Nthreads);
	 omp_set_num_threads(run.Nthreads);
	#endif

	Ngrid=run.Ngrid;
	Ngrid3=Ngrid*Ngrid*Ngrid; //total number of grid points.
	Header=load_snapshot_header(0);
	boxsize=Header.BoxSize;
	gridsize=boxsize/Ngrid;
	OmM=Header.Omega0;
	OmL=Header.OmegaLambda;
	Header=load_galaxies_header(0);
	NumGal=Header.npartTotal[1];
	printf("	Box: %lf	GravGrid: %i	Gridsize: %lf	#Particles: %i OmM: %f	OmL: %f\n",boxsize,Ngrid,gridsize,NumGal,Header.Omega0,Header.OmegaLambda);
	printf("	Observer centred on:	%lf	%lf	%lf\n",run.ObsPos[0],run.ObsPos[1],run.ObsPos[2]);

	if (!(Head=calloc(run.NSnap,sizeof(struct header_snapshot)))) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
	if (!(rangemin=calloc(run.NSnap,sizeof(double)))) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
	if (!(rangemax=calloc(run.NSnap,sizeof(double)))) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
	if (!(PotHold=malloc(Ngrid3*sizeof(float)))) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
	if (!(Hubarray=malloc(run.NSnap*sizeof(double)))) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
	if (!(Pot=calloc(Ngrid3,sizeof(void*)))) {printf("Failed to allocate memory for the grid.\n"); exit(0);}

	NMaxOutput=abs(run.OutBuf*NumGal*4.*CONSPI/3.*(run.maxdist*run.maxdist*run.maxdist-run.mindist*run.mindist*run.mindist)/(boxsize*boxsize*boxsize));
	if (NMaxOutput > NumGal) NMaxOutput=NumGal; //Estimate how many particles will be in the output array.
	printf("	Reserve space for %i output particles.\n",NMaxOutput);

	if (!(AllConePos=calloc(NMaxOutput,sizeof(struct light_cone_object)))) {printf("Failed to allocate memory for the grid.\n"); exit(0);}

	#ifdef MEASURETIME
	 PartTime[0]-=time(NULL);
	#endif
	//Setting up the potential etc for the integration.
	printf("\n\n.....................................................\nSetting-up for the light cone integration\n");
	fflush(stdout);
	for (i=0;i<run.NSnap;i++){//Internally we want to have the latest snapshot as the first in the array.
		printf("Processing snapshot %i\n",run.NSnap-1-i);
		if (!(Pot[i]=malloc(run.NSnap*sizeof(struct potential_container)))) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
		Head[i]=load_snapshot_header(run.NSnap-1-i);
		printf("	z=%f\n",Head[i].redshift);
		fflush(stdout);
		if (fabs(Head[i].redshift) < 1e-8) rangemax[i]=0; else rangemax[i]=calc_comoving_distance(Head[i].redshift,Head[i].redshift/5000,OmM,OmL);
		Hubarray[i]=Hub(Head[i].time,OmM,OmL)*Head[i].time; //Conformal time.
		if (i!=0) rangemin[i]=rangemax[i-1]; else rangemin[i]=rangemax[i];

		sprintf(gravpath,"%s/GravPot%03d_FullSize_Ng%1d_cic.dat",run.gravpath,run.NSnap-1-i,Ngrid);
		get_potential(gravpath,&PotHold,Ngrid,Head[i],path); //Read the potential from a file or compute it.

		fill_potential_container(PotHold,Ngrid,boxsize,Pot[i]); //Computes all spacial derivatives.
	}
	compute_potential_timeder(Ngrid3,run.NSnap,rangemax,Head,Pot); //Computes all time derivatives.

	free(PotHold);//free unused memory.
	#ifdef MEASURETIME
	 PartTime[1]-=time(NULL);
	 PartTime[0]-=PartTime[1];
	#endif

	//omp_set_num_threads(1);//Debugging mode only.
	d=0;
	printf("\n\n.....................................................\n");
	printf("Calculating displacement ... 000");
	fflush(stdout);
	#pragma omp parallel private(a,b,c,j,ray,AllP,RealStepPart,bStart,bStop,ConePos)
	{
	#pragma omp critical
	{ //each thread allocates some working memory.
	 RealStepPart=run.NPartInMem;
	 if (!(ray=malloc(run.NSnap*sizeof(struct light_rays_path)))) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
	 if (!(AllP=malloc(run.NSnap*sizeof(void*)))) {printf("Failed to allocate memory for the particles.\n"); exit(0);}
	}
	#pragma omp for schedule(dynamic,1)
	for (c=0;c<NumGal;c+=run.NPartInMem){
		//Reading run.NPartInMem particles at all available times!
		for (j=0;j<run.NSnap;j++){//Internally we want to have the latest snapshot as the first in the array.
			RealStepPart=run.NPartInMem; //TODO: Push this to load_galaxies_partial and output the particle array.
			if (c+run.NPartInMem > NumGal) RealStepPart=NumGal-c;
			#pragma omp critical (Reading)
			{
			 while (load_galaxies_partial(run.NSnap-1-j,1,c,RealStepPart,&AllP[j])!=0) {
				sleep(1);
			 }
			}
		}
		for (a=0;a<RealStepPart;a++) {

			for (b=0;b<run.NSnap;++b) {
				ray[b].D3shift[0]=0; ray[b].D3shift[1]=0; ray[b].D3shift[2]=0; //We will explicitly ask for this to ensure the shift has been calculated.
				initialize_ray(&ray[b],AllP[b][a],Head[b].BoxSize,rangemax[b],Hubarray[b]);
			}

			//for (b=0;b<run.NSnap;++b) printf("	%f	%f	%f	%f	%f	%f	|	%f	%f	%f	|	%f	%f	%f	%f\n",ray[b].start[0],ray[b].start[1],ray[b].start[2],ray[b].end[0],ray[b].end[1],ray[b].end[2],PeriDist(ray[b].end[0]-ray[b].start[0],boxsize),PeriDist(ray[b].end[1]-ray[b].start[1],boxsize),PeriDist(ray[b].end[2]-ray[b].start[2],boxsize),ray[b].dist,ray[b].vrHub,rangemax[b],ray[b].TimeDistance);
			//estimates which line-of-sights need to be integrated, and if target is inside the range
			if (get_integration_range(ray,rangemin,rangemax,gridsize,&bStart,&bStop)) continue; 

			for (b=bStart;b<bStop;++b){
				ray[b].end[0]+=ray[b].vrHub*ray[b].normvec[0];//Shift position by the particles radial velocity.
				ray[b].end[1]+=ray[b].vrHub*ray[b].normvec[1];
				ray[b].end[2]+=ray[b].vrHub*ray[b].normvec[2];
				ray[b].TimeDistance-=ray[b].vrHub;
				update_ray(&ray[b]);

				ray[b].pot=get_potential_value(ray[b].end[0],ray[b].end[1],ray[b].end[2],gridsize,Ngrid,Pot,b);

				find_cells_along_ray(&ray[b],gridsize,Ngrid);
				if (ray[b].CellsAlongPath == 0) continue;

				integrate_path(&ray[b],Hubarray[b],rangemin,rangemax,Pot);

				free(ray[b].PathIntersec); //Free memory allocated by find_cells_along_ray
				free(ray[b].PathCell);

				ray[b].end[0]-=ray[b].DRshift*ray[b].normvec[0]+ray[b].D3shift[0]; //Shift position by the remaining shifts.
				ray[b].end[1]-=ray[b].DRshift*ray[b].normvec[1]+ray[b].D3shift[0];
				ray[b].end[2]-=ray[b].DRshift*ray[b].normvec[2]+ray[b].D3shift[0];
				ray[b].TimeDistance-=ray[b].DTshift;
				update_ray(&ray[b]);
			}

			if (get_integration_range(ray,rangemin,rangemax,gridsize,&bStart,&bStop)) continue; //Required range might change due to the shift.
			
			for (b=bStart;b<bStop;++b){
				if (ray[b].D3shift[0]!=0.0) continue; //Shift has already been computed.
				ray[b].end[0]+=ray[b].vrHub*ray[b].normvec[0];//Shift position by the particles radial velocity.
				ray[b].end[1]+=ray[b].vrHub*ray[b].normvec[1];
				ray[b].end[2]+=ray[b].vrHub*ray[b].normvec[2];
				ray[b].TimeDistance-=ray[b].vrHub;
				update_ray(&ray[b]);

				ray[b].pot=get_potential_value(ray[b].end[0],ray[b].end[1],ray[b].end[2],gridsize,Ngrid,Pot,b);

				find_cells_along_ray(&ray[b],gridsize,Ngrid);
				if (ray[b].CellsAlongPath == 0) continue;

				integrate_path(&ray[b],Hubarray[b],rangemin,rangemax,Pot);
				
				free(ray[b].PathIntersec); //Free memory allocated by find_cells_along_ray
				free(ray[b].PathCell);

				ray[b].end[0]-=ray[b].DRshift*ray[b].normvec[0]+ray[b].D3shift[0]; //Shift position by the remaining shifts.
				ray[b].end[1]-=ray[b].DRshift*ray[b].normvec[1]+ray[b].D3shift[0];
				ray[b].end[2]-=ray[b].DRshift*ray[b].normvec[2]+ray[b].D3shift[0];
				ray[b].TimeDistance-=ray[b].DTshift;
				update_ray(&ray[b]);
			}

			if (get_all_light_cone_intersection(ray,Head[0],&ConePos)==0) {
				#pragma omp critical (NextPart)
				{ //Make sure only one thead writes to the outut array at a time.
				if (NextPart==NMaxOutput) {
					printf("\nNot enough space in output array.\nIncrease MemBufFactor and restart.\n");
					exit(0);
				}
				//ConePos.id=c+a;
				AllConePos[NextPart]=ConePos;
				++NextPart;
				if (c+a==989761) printf("%i\n",NextPart);
				}
			}
		}
		#pragma omp atomic
		d+=RealStepPart;
		printf("\b\b\b%03d",abs(d/((double)NumGal)*100));
		fflush(stdout);
	for (j=0;j<run.NSnap;++j) free(AllP[j]);
	}
	free(ray);
	free(AllP);
	}
	printf("done.\nNumber of particles on the light cone: %i\n",NextPart);
	fflush(stdout);
	#ifdef MEASURETIME
	 PartTime[1]+=time(NULL);
	#endif

	LCsize=2*run.maxdist;
	sprintf(path,"%s/%s_realspace.dat",run.outpath,run.outname);
	printf("Writing to %s ... ",path);
	fflush(stdout);
	if ((out=fopen(path,"w"))) {
		fwrite(&LCsize,sizeof(float),1,out);
		fwrite(&NextPart,sizeof(int),1,out);
		for (i=0;i<NextPart;++i) fwrite(&AllConePos[i].SimPos[0],sizeof(float),3,out);
		fclose(out);
		printf("done.\n");
	} else printf("Can not write to path!\n");
	sprintf(path,"%s/%s_vRSD.dat",run.outpath,run.outname);
	printf("Writing to %s ... ",path);
	fflush(stdout);
	if ((out=fopen(path,"w"))) {
		fwrite(&LCsize,sizeof(float),1,out);
		fwrite(&NextPart,sizeof(int),1,out);
		for (i=0;i<NextPart;++i) fwrite(&AllConePos[i].VelPos[0],sizeof(float),3,out);
		fclose(out);
		printf("done.\n");
	} else printf("Can not write to path!\n");
	sprintf(path,"%s/%s_GRRSD.dat",run.outpath,run.outname);
	printf("Writing to %s ... ",path);
	fflush(stdout);
	if ((out=fopen(path,"w"))) {
		fwrite(&LCsize,sizeof(float),1,out);
		fwrite(&NextPart,sizeof(int),1,out);
		for (i=0;i<NextPart;++i) fwrite(&AllConePos[i].RelPos[0],sizeof(float),3,out);
		fclose(out);
		printf("done.\n");
	} else printf("Can not write to path!\n");
	sprintf(path,"%s/%s_magni.dat",run.outpath,run.outname);
	printf("Writing to %s ... ",path);
	fflush(stdout);
	if ((out=fopen(path,"w"))) {
		fwrite(&LCsize,sizeof(float),1,out);
		fwrite(&NextPart,sizeof(int),1,out);
		for (i=0;i<NextPart;++i) fwrite(&AllConePos[i].Magni[0],sizeof(float),1,out);
		for (i=0;i<NextPart;++i) fwrite(&AllConePos[i].Magni[1],sizeof(float),1,out);
		for (i=0;i<NextPart;++i) fwrite(&AllConePos[i].Magni[2],sizeof(float),1,out);
		fclose(out);
		printf("done.\n");
	} else printf("Can not write to path!\n");

	#ifdef MEASURETIME
	 TotalTime+=time(NULL);
	 sprintf(path,"%s/%s_runtime.txt",run.outpath,run.outname);
	 printf("Writing to %s ... ",path);
	 fflush(stdout);
	 if ((out=fopen(path,"w"))) {
		fprintf(out,"Total wall-clock time in sec: %.2f\n",(double)TotalTime);
		fprintf(out," - Potential calculation: %.2f\n",(double)PartTime[0]);
		fprintf(out,"    - including I/O: %.2f\n",(double)IOTime[0]);
		fprintf(out," - Main loop: %.2f\n",(double)PartTime[1]);
		fprintf(out,"    - including I/O (per thread): %.2f\n",(double)IOTime[1]/run.Nthreads);
		fclose(out);
		printf("done.\n");
	 } else printf("Can not write to path!\n");
	#endif

	return 0;
}

