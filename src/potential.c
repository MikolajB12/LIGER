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

//  This file contains the routines which compute the potential and its derivatives.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>   
#include "./vardef.h"
#include "./load_snap.h"
#include "./potential.h"

void fill_potential_container(float *PotInput,int Ngrid,double boxsize,struct potential_container *PotCont){ //Computes all spacial derivatives.
double dum,kscale;
int a,b,c,i,koor,k1,k2,k3,Ngrid3=Ngrid*Ngrid*Ngrid;
double *pot;
fftw_complex *pot_fft,*der1_fft,*der2_fft,*der3_fft;
fftw_plan PotFFT,PotDer1FFT,PotDer2FFT,PotDer3FFT;

printf("	Computing fft ...");
fflush(stdout);

if (!(pot=(double*) fftw_malloc(Ngrid3*sizeof(double)))) {printf("\nAllocation of memory failed.\n"); exit(0);}
if (!(pot_fft=(fftw_complex*) fftw_malloc(Ngrid*Ngrid*(Ngrid/2+1)*sizeof(fftw_complex)))) {printf("\nAllocation of memory failed.\n"); exit(0);}
if (!(der1_fft=(fftw_complex*) fftw_malloc(Ngrid*Ngrid*(Ngrid/2+1)*sizeof(fftw_complex)))) {printf("\nAllocation of memory failed.\n"); exit(0);}
if (!(der2_fft=(fftw_complex*) fftw_malloc(Ngrid*Ngrid*(Ngrid/2+1)*sizeof(fftw_complex)))) {printf("\nAllocation of memory failed.\n"); exit(0);}
if (!(der3_fft=(fftw_complex*) fftw_malloc(Ngrid*Ngrid*(Ngrid/2+1)*sizeof(fftw_complex)))) {printf("\nAllocation of memory failed.\n"); exit(0);}

#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) {
	pot[i]=PotInput[i]; //Need to change to double for usage of fftw3.
	PotCont[i].val=PotInput[i];
}

PotFFT=fftw_plan_dft_r2c_3d(Ngrid, Ngrid, Ngrid, pot, pot_fft, FFTW_ESTIMATE);
PotDer1FFT=fftw_plan_dft_c2r_3d(Ngrid, Ngrid, Ngrid, der1_fft, pot, FFTW_ESTIMATE);
PotDer2FFT=fftw_plan_dft_c2r_3d(Ngrid, Ngrid, Ngrid, der2_fft, pot, FFTW_ESTIMATE);
PotDer3FFT=fftw_plan_dft_c2r_3d(Ngrid, Ngrid, Ngrid, der3_fft, pot, FFTW_ESTIMATE);
fftw_execute(PotFFT); /* Fourier transform */

kscale=2.0*CONSPI/(boxsize*Ngrid3);
printf("done.\n	Computing 1st derivative ...");
fflush(stdout);
#pragma omp parallel for schedule(guided) private(k1,k2,k3,dum,a,b,c,koor) firstprivate(kscale)
for (a=0;a<Ngrid;a++){
	for (b=0;b<Ngrid;b++){
		for (c=0;c<(Ngrid/2+1);c++){ 
			k1=(a>(Ngrid-1)/2 ? (-Ngrid+a):a);
			k2=(b>(Ngrid-1)/2 ? (-Ngrid+b):b);
			k3=c;
			koor=c+(Ngrid/2+1)*(b+Ngrid*a);

			if (k1!=0) dum=k1*kscale; else dum=kscale;
			der1_fft[koor][0]=-dum*pot_fft[koor][1];
			der1_fft[koor][1]=dum*pot_fft[koor][0];

			if (k2!=0) dum=k2*kscale; else dum=kscale;
			der2_fft[koor][0]=-dum*pot_fft[koor][1];
			der2_fft[koor][1]=dum*pot_fft[koor][0];

			if (k3!=0) dum=k3*kscale; else dum=kscale;
			der3_fft[koor][0]=-dum*pot_fft[koor][1];
			der3_fft[koor][1]=dum*pot_fft[koor][0];
		}
	}
}

fftw_execute(PotDer1FFT); /* Fourier transform back */
#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) PotCont[i].firD[0]=pot[i];
fftw_execute(PotDer2FFT); /* Fourier transform back */
#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) PotCont[i].firD[1]=pot[i];
fftw_execute(PotDer3FFT); /* Fourier transform back */
#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) PotCont[i].firD[2]=pot[i];


kscale=4.0*CONSPI*CONSPI/(boxsize*boxsize*Ngrid3);
printf("done.\n	Computing 2nd derivative ...");
fflush(stdout);
#pragma omp parallel for schedule(guided) private(k1,k2,k3,dum,a,b,c,koor) firstprivate(kscale)
for (a=0;a<Ngrid;a++){
	for (b=0;b<Ngrid;b++){
		for (c=0;c<(Ngrid/2+1);c++){ 
			k1=(a>(Ngrid-1)/2 ? (-Ngrid+a):a);
			k2=(b>(Ngrid-1)/2 ? (-Ngrid+b):b);
			k3=c;
			koor=c+(Ngrid/2+1)*(b+Ngrid*a);

			//dum=pow(k1,2.);//+pow(k2,2.)+pow(k3,2.);
			//if (dum==0) dum=-1; else dum*=(-1);
			//dum*=2.0*CONSPI/boxsize*2.0*CONSPI/boxsize;
			//if ((a==Ngrid-3)&&(b==Ngrid-1)&&(c==1)) printf("%i %i %i %lf\n%lf ",a,b,c,dum,pot_fft[koor][0]);
			if (k1!=0) dum=-k1*k1*kscale; else dum=-kscale;
			der1_fft[koor][0]=dum*pot_fft[koor][0];
			der1_fft[koor][1]=dum*pot_fft[koor][1];

			if (k2!=0) dum=-k2*k2*kscale; else dum=-kscale;
			der2_fft[koor][0]=dum*pot_fft[koor][0];
			der2_fft[koor][1]=dum*pot_fft[koor][1];

			if (k3!=0) dum=-k3*k3*kscale; else dum=-kscale;
			der3_fft[koor][0]=dum*pot_fft[koor][0];
			der3_fft[koor][1]=dum*pot_fft[koor][1];
		}
	}
}

fftw_execute(PotDer1FFT); /* Fourier transform back */
#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) PotCont[i].secD[0]=pot[i];
fftw_execute(PotDer2FFT); /* Fourier transform back */
#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) PotCont[i].secD[1]=pot[i];
fftw_execute(PotDer3FFT); /* Fourier transform back */
#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) PotCont[i].secD[2]=pot[i];

#pragma omp parallel for schedule(guided) private(k1,k2,k3,dum,a,b,c,koor) firstprivate(kscale)
for (a=0;a<Ngrid;a++){
	for (b=0;b<Ngrid;b++){
		for (c=0;c<(Ngrid/2+1);c++){
			k1=(a>(Ngrid-1)/2 ? (-Ngrid+a):a);
			k2=(b>(Ngrid-1)/2 ? (-Ngrid+b):b);
			k3=c;
			koor=c+(Ngrid/2+1)*(b+Ngrid*a);

			//dum=k1*k2;
			//if (dum==0) dum=-1; else dum*=(-1);
			//dum*=2.0*CONSPI/boxsize*2.0*CONSPI/boxsize;
			if (k1*k2!=0) dum=-k1*k2*kscale; else dum=-kscale;
			der1_fft[koor][0]=dum*pot_fft[koor][0];
			der1_fft[koor][1]=dum*pot_fft[koor][1];

			if (k1*k3!=0) dum=-k1*k3*kscale; else dum=-kscale;
			der2_fft[koor][0]=dum*pot_fft[koor][0];
			der2_fft[koor][1]=dum*pot_fft[koor][1];

			if (k2*k3!=0) dum=-k2*k3*kscale; else dum=-kscale;
			der3_fft[koor][0]=dum*pot_fft[koor][0];
			der3_fft[koor][1]=dum*pot_fft[koor][1];
		}
	}
}

fftw_free(pot_fft);

fftw_execute(PotDer1FFT); /* Fourier transform back */
#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) PotCont[i].croD[0]=pot[i];
fftw_execute(PotDer2FFT); /* Fourier transform back */
#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) PotCont[i].croD[1]=pot[i];
fftw_execute(PotDer3FFT); /* Fourier transform back */
#pragma omp parallel for schedule(guided)
for (i=0;i<Ngrid3;i++) PotCont[i].croD[2]=pot[i];

fftw_free(der1_fft);
fftw_free(der2_fft);
fftw_free(der3_fft);

fftw_destroy_plan(PotFFT);
fftw_destroy_plan(PotDer1FFT);
fftw_destroy_plan(PotDer2FFT);
fftw_destroy_plan(PotDer3FFT);

printf("done.\n");
fflush(stdout);
}


/*Here the overdensity field is calculated from the density field*/
int overdens(double *Alldens,int Ngrid){
	long a;
	double sum=0;
	printf("	Rescale density to contrast ... ");
#pragma omp parallel for reduction(+:sum)
	for (a=0;a<Ngrid*Ngrid*Ngrid;a++) sum+=Alldens[a];
	sum/=(Ngrid*Ngrid*Ngrid);
	//printf("Total density: %lf 10^10M_sol\n",sum);
	//printf("Mean density rescaled: %le -> ",sum);
#pragma omp parallel for schedule(guided)
	for (a=0;a<Ngrid*Ngrid*Ngrid;a++) Alldens[a]=(Alldens[a]-sum)/sum;
#pragma omp parallel for reduction(+:sum)
	for (a=0;a<Ngrid*Ngrid*Ngrid;a++) sum+=Alldens[a];
	sum/=(Ngrid*Ngrid*Ngrid);
	printf("done.\n");
return 0;
}



// Here the density field is calculated on a grid Ngrid^3 
double *density_field(struct particle_data *P,struct header_snapshot header1,unsigned long Ngrid,int flag){
long i,xkoor,ykoor,zkoor;
int a,b,c,NumPart;
double gridsize,fak;
double *Alldens;
printf("	Computing density ");
fflush(stdout);
if (!(Alldens=calloc(Ngrid*Ngrid*Ngrid,sizeof(double)))){
	printf("Failed to allocate memory for the grid.\n");
	exit(0);
}
NumPart=header1.npartTotal[1];
gridsize=header1.BoxSize/Ngrid;
if (flag>0){
	printf("with CIC ... ");
	fflush(stdout);
	#pragma omp parallel for schedule(guided) private(a,b,c,xkoor,ykoor,zkoor)
	for (i=0;i<NumPart;i++){
		xkoor=abs(P[i].Pos[0]/gridsize);
		ykoor=abs(P[i].Pos[1]/gridsize);
		zkoor=abs(P[i].Pos[2]/gridsize);
		for (a=-1;a<2;a++){
			for (b=-1;b<2;b++){
				for (c=-1;c<2;c++){
					#pragma omp atomic
					#ifdef NOUNIGRID
					 Alldens[shift_koor(xkoor+a,ykoor+b,zkoor+c,Ngrid)]+=P[i].Mass*weight(P[i],(xkoor+0.5+a)*gridsize,(ykoor+0.5+b)*gridsize,(zkoor+0.5+c)*gridsize,gridsize,1);
					#else
					 Alldens[shift_koor(xkoor+a,ykoor+b,zkoor+c,Ngrid)]+=weight(P[i],(xkoor+0.5+a)*gridsize,(ykoor+0.5+b)*gridsize,(zkoor+0.5+c)*gridsize,gridsize,1);
					#endif
				}
			}
		}
	}
} else {
	printf("with NGB ... ");
	fflush(stdout);
	#pragma omp parallel for schedule(guided) private(a,b,c,xkoor,ykoor,zkoor)
	for (i=0;i<NumPart;i++){
		xkoor=abs(P[i].Pos[0]/gridsize);
		ykoor=abs(P[i].Pos[1]/gridsize);
		zkoor=abs(P[i].Pos[2]/gridsize);
		#pragma omp atomic
		#ifdef NOUNIGRID
		 Alldens[shift_koor(xkoor,ykoor,zkoor,Ngrid)]+=P[i].Mass;
		#else
		 Alldens[shift_koor(xkoor,ykoor,zkoor,Ngrid)]++;
		#endif

	}
}
#ifdef NOUNIGRID
 fak=pow(Ngrid/header1.BoxSize,3.);
#else
 fak=header1.mass[1]*pow(Ngrid/header1.BoxSize,3.);
#endif
for (a=0;a<Ngrid*Ngrid*Ngrid;a++) Alldens[a]*=fak;

printf("done.\n");
fflush(stdout);

return Alldens;
}

/* This function computes the particle weights for the CIC grid assignment */
double weight(struct particle_data part,double gridx,double gridy,double gridz,double gridsize,int flag){
	double dum=1.,dist;
	double grid[3]={gridx,gridy,gridz};
	int w;
	if (flag == 1) for (w=0;w<3;w++){
		dist=(1.-fabs(part.Pos[w]-grid[w])/gridsize);
		if (dist<0) return 0;
		dum*=dist;
		}
	else if (flag == 2) for (w=0;w<3;w++){
				if (fabs(part.Pos[w]-grid[w]) < gridsize/2.) dum*=(3./4.-(part.Pos[w]-grid[w])*(part.Pos[w]-grid[w])/gridsize/gridsize);
				else if (fabs(part.Pos[w]-grid[w]) < 3.*gridsize/2.) dum*=(3./2.-fabs(part.Pos[w]-grid[w])/gridsize)*(3./2.-fabs(part.Pos[w]-grid[w])/gridsize)/2.;
				else return 0;
				}
	if (dum > 0) return dum;
	else return 0;
}

//This function computes the potential from the density (contrast).
float *grav_potential(double **dens,int Ngrid,struct header_snapshot header1){
double dum,k1,k2,k3,kN,gridsize;
int a,b,c,koor;
double *pot;
float *result;
fftw_complex *dens_fft;
fftw_plan poisson;

printf("	Computing grav. potential:\n	FFT of Density field ... ");
fflush(stdout);
if (!(dens_fft=(fftw_complex*) fftw_malloc(Ngrid*Ngrid*(Ngrid/2+1)*sizeof(fftw_complex)))){
	printf("\nAllocation of memory failed.\n");
	exit(0);
}

gridsize=header1.BoxSize/Ngrid;
kN=CONSPI/gridsize; //Nyquist freq
poisson=fftw_plan_dft_r2c_3d(Ngrid, Ngrid, Ngrid, dens[0], dens_fft, FFTW_ESTIMATE);
fftw_execute(poisson); /* Fourier transform */
printf("done.\n	Solving Poisson ... ");
fflush(stdout);
#pragma omp parallel for schedule(dynamic,1) private(k1,k2,k3,dum,a,b,c,koor)
for (a=0;a<Ngrid;a++){
	for (b=0;b<Ngrid;b++){
		for (c=0;c<(Ngrid/2+1);c++){ /* Solving poisson's equation in fourier space */
			if (a>(Ngrid-1)/2) k1=((double)(-Ngrid+a))/(Ngrid*gridsize)*2.0*CONSPI; else k1=a/(Ngrid*gridsize)*2.0*CONSPI;
			if (b>(Ngrid-1)/2) k2=((double)(-Ngrid+b))/(Ngrid*gridsize)*2.0*CONSPI; else k2=b/(Ngrid*gridsize)*2.0*CONSPI;
			k3=c/(Ngrid*gridsize)*2.0*CONSPI;
			//if (a>(Ngrid-1)/2) k1=(-Ngrid+a)*2.0/(2*Ngrid+1)*2.0*CONSPI; else k1=a*2.0/(2*Ngrid+1)*2.0*CONSPI;
			//if (b>(Ngrid-1)/2) k2=(-Ngrid+b)*2.0/(2*Ngrid+1)*2.0*CONSPI; else k2=b*2.0/(2*Ngrid+1)*2.0*CONSPI;
			//k3=c*2.0/(2*Ngrid+1)*2.0*CONSPI;
			dum=(k1*k1+k2*k2+k3*k3);//*pow(sin(0.5*k1/kN)*sin(0.5*k2/kN)*sin(0.5*k3/kN)/(0.125*k1*k2*k3/kN/kN/kN),2.); //Correcting for the CIC
			if (k1!=0) dum*=pow(sin(CONSPI*0.5*k1/kN)/(CONSPI*0.5*k1/kN),2.);
			if (k2!=0) dum*=pow(sin(CONSPI*0.5*k2/kN)/(CONSPI*0.5*k2/kN),2.);
			if (k3!=0) dum*=pow(sin(CONSPI*0.5*k3/kN)/(CONSPI*0.5*k3/kN),2.);
			//if (dum==0) {printf("Zero window function!\n%i	%i	%i\n%e	%e	%e\n%e	%e	%e",a,b,c,k1,k2,k3,sin(0.5*k1/kN),sin(0.5*k2/kN),sin(0.5*k3/kN)); exit(0);}
				koor=c+(Ngrid/2+1)*(b+Ngrid*a);
			if (dum==0) {
				dens_fft[koor][0]=0.0;
				dens_fft[koor][1]=0.0;
			} else {
				dens_fft[koor][0]/=-dum;
				dens_fft[koor][1]/=-dum;
			}
			if ((isfinite(dens_fft[koor][0])==0)||(isfinite(dens_fft[koor][1])==0)) {printf("Non-finite entry!\n%e	%e	%e\n%i	%i	%i\n%e	%e	%e\n%e	%e	%e",dens_fft[koor][0],dens_fft[koor][1],dum,a,b,c,k1,k2,k3,sin(0.5*k1/kN),sin(0.5*k2/kN),sin(0.5*k3/kN)); exit(0);}
			}
		}
	}
printf("done.\n	FFT to real space ... ");
fflush(stdout);
if (!(pot=malloc(Ngrid*Ngrid*Ngrid*sizeof(double)))){
	printf("\nAllocation of memory failed.\n");
	exit(0);
}
poisson=fftw_plan_dft_c2r_3d(Ngrid, Ngrid, Ngrid, dens_fft, pot, FFTW_ESTIMATE);
fftw_execute(poisson); /* Fourier transform back */
fftw_free(dens_fft);
if (!(result=malloc(Ngrid*Ngrid*Ngrid*sizeof(float)))){
	printf("\nAllocation of memory failed.\n");
	exit(0);
}
dum=1./(Ngrid*Ngrid*Ngrid);
#pragma omp parallel for schedule(guided)
for (a=0;a<Ngrid*Ngrid*Ngrid;a++) result[a]=pot[a]*dum*1.5*100.0*100.0*header1.Omega0/header1.time;// Normalizing: this now has km^2/s^2 as unit.
printf("done.\n");
fflush(stdout);
return result;
}

double get_potential_value(double posx,double posy,double posz,double gridsize,int Ngrid,struct potential_container *Pot){
int x,y,z;
x=abs(posx/gridsize);
y=abs(posy/gridsize);
z=abs(posz/gridsize);
return Pot[shift_koor(x,y,z,Ngrid)].val;
}

void get_potential(char *fname,float **pot,int Ngrid,struct header_snapshot Head,char *path){
FILE *fileptr;
int Ncheck;
float *PotHold=pot[0];
double *dens;
struct particle_data *Pptr;

if ((fileptr=fopen(fname,"r"))) {
	#ifdef MEASURETIME
	 IOTime[0]-=time(NULL);
	#endif
	printf("	Graviational potential: \n	Reading from %s ... ",fname);
	fflush(stdout);
	fread(&Ncheck,sizeof(int),1,fileptr);
	if (Ncheck!=Ngrid) {printf("File contains unexpected grid resolution!\nFile	Desired\n%i	%i\nRemove file and restart.",Ncheck,Ngrid); exit(0);}
	fread(PotHold,sizeof(float),Ncheck*Ncheck*Ncheck,fileptr);
	fclose(fileptr);
	printf("done.\n");
	fflush(stdout);
	#ifdef MEASURETIME
	 IOTime[0]+=time(NULL);
	#endif
} else {
	//Here we then compute the potential from the particle distribution.
	printf("	Graviational potential: \n	File not found: %s -> Calculating:\n",fname);
	fflush(stdout);
	load_snapshot(path,1,&Pptr);
	dens=density_field(Pptr,Head,Ngrid,1);
	free(Pptr);
	overdens(dens,Ngrid);
	PotHold=grav_potential(&dens,Ngrid,Head);
	if (run.ifOutputPotential) {
		if ((fileptr=fopen(fname,"w"))) {
			printf("	Saving potantial ... ");
			fflush(stdout);
			fwrite(&Ngrid,sizeof(int),1,fileptr);
			fwrite(PotHold,sizeof(float),Ngrid*Ngrid*Ngrid,fileptr);
			fclose(fileptr);
			printf("done.\n");
			fflush(stdout);
		} else {
			printf("	Cannot open file to save potential!\n");
			fflush(stdout);
		}
	}
}
}

void compute_potential_timeder(int Ngrid3,int Nsnap,double *rangemax,struct header_snapshot *Head,struct potential_container **PotCont){
int i,koor;
double weight,dist;

printf("Computing first time derivative ... ");
fflush(stdout);
for (i=1;i<Nsnap-1;i++){
	dist=(rangemax[i+1]-rangemax[i])/(rangemax[i]-rangemax[i-1]);
	weight=(rangemax[i]-rangemax[i-1])+(rangemax[i+1]-rangemax[i])*(rangemax[i+1]-rangemax[i])/(rangemax[i]-rangemax[i-1]);
	dist=dist*dist;
	for (koor=0;koor<Ngrid3;++koor) PotCont[i][koor].timeD=(PotCont[i+1][koor].val-PotCont[i-1][koor].val*dist)/weight*Head[i].time;
}
weight=1.0/(rangemax[1]-rangemax[0]);
dist=1.0/(rangemax[Nsnap-1]-rangemax[Nsnap-2]);
#pragma omp parallel for schedule(guided) private(koor)
for (koor=0;koor<Ngrid3;koor++){
	PotCont[0][koor].timeD=(PotCont[1][koor].val-PotCont[0][koor].val)*weight*Head[0].time;
	PotCont[Nsnap-1][koor].timeD=(PotCont[Nsnap-1][koor].val-PotCont[Nsnap-2][koor].val)*dist*Head[Nsnap-1].time;
}
printf("done.\n");
fflush(stdout);
}


int shift_koor(int a,int b,int c,int Ngrid){
int dum,mem[3];
mem[0]=a;
mem[1]=b;
mem[2]=c;
if (a>=Ngrid){dum=a/Ngrid; a=a-Ngrid*dum;} else if (a<0){dum=abs((-a-1)/Ngrid)+1; a=Ngrid*dum+a;}
if (b>=Ngrid){dum=b/Ngrid; b=b-Ngrid*dum;} else if (b<0){dum=abs((-b-1)/Ngrid)+1; b=Ngrid*dum+b;}
if (c>=Ngrid){dum=c/Ngrid; c=c-Ngrid*dum;} else if (c<0){dum=abs((-c-1)/Ngrid)+1; c=Ngrid*dum+c;}
if ((a>=Ngrid)||(b>=Ngrid)||(c>=Ngrid)||(a<0)||(b<0)||(c<0)) {printf("ERROR:SHIFT_KOOR:%i (MAX=%i/a=%i-%i/b=%i-%i/c=%i-%i)\n",c+Ngrid*(b+Ngrid*a),Ngrid*Ngrid*Ngrid,mem[0],a,mem[1],b,mem[2],c); exit(0);}
return c+Ngrid*(b+Ngrid*a);
}

double get_distance_square_periodic(double x1,double x2,double x3,double y1,double y2,double y3,double boxsize){
double dist;
	dist=0;
	if (fabs(x1-y1)>boxsize/2.) dist+=pow(boxsize-fabs(x1-y1),2.);
	else dist+=pow(x1-y1,2.);
	if (fabs(x2-y2)>boxsize/2.) dist+=pow(boxsize-fabs(x2-y2),2.);
	else dist+=pow(x2-y2,2.);
	if (fabs(x3-y3)>boxsize/2.) dist+=pow(boxsize-fabs(x3-y3),2.);
	else dist+=pow(x3-y3,2.);
	return dist;
}

double PeriDist(double dist,double box){
//int dum;
//if (dist>=box){dum=dist/box; return dist-box*dum;} else if (dist<0){dum=abs((-dist-1)/box)+1; return box*dum+dist;}
if (dist>box/2.) return dist-box; else if (dist<-box/2.) return dist+box; else return dist;
}

