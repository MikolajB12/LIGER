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

//  This file contains the routines to perform the integration along the line of sight.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>  
#include "./vardef.h"

double Hub(double a,double Om,double OL){
return 100.0*sqrt(Om/a/a/a+(1.0-Om-OL)/a/a+OL);
//hub=(0.0001022729969);// in units of 1/Myr
}

double calc_comoving_distance(double red,double dz,double Om,double OL){
int i=0;
double z,result=0.0;
result=pow(Om*pow(red+1,3)+(red+1)*(red+1)*(1-Om-OL)+OL,-0.5)/2.0+0.5;
for (z=dz;z<red-dz*0.1;z+=dz,i++) result+=pow(Om*pow(z+1,3)+(z+1)*(z+1)*(1-Om-OL)+OL,-0.5);
return result*dz*2997.92458; // Mpc/h
}

//interpolates four points to the desired location.
double simulation_time_interpol(double fieldcloser,double fieldclose,double fieldfar,double fieldfarer,double distcloser,double distclose,double distfar,double distfarer,double dist){
double pg,g,p,c,d,dfdxc,dfdxf;
g=distfar-distclose;
p=(dist-distclose);
if ((p < 0)||(p > g)) {
	#pragma omp critical
	{
	printf("Error, out of range interpolation!\nWanted:	%le	(%le	%le)\nBoundary:	%le	%le	%le	%le\nField:	%le	%le	%le	%le\n",dist,p,g,distcloser,distclose,distfar,distfarer,fieldcloser,fieldclose,fieldfar,fieldfarer);
	fflush(stdout);
	exit(0);
	}
}
if (fieldfarer == fieldfar) dfdxf=(fieldfar-fieldclose)/g;
else {
	c=g/(distfarer-distfar);
	dfdxf=(fieldclose-fieldfar-(fieldfarer-fieldfar)*c*c)/(-g*(1.0+c));
}
if (fieldcloser == fieldclose) dfdxc=(fieldfar-fieldclose)/g;
else {
	c=(distcloser-distclose)/g;
	dfdxc=(fieldcloser-fieldclose-(fieldfar-fieldclose)*c*c)/((distcloser-distclose)*(1.0-c));
}
pg=p/g;
c=(3*(fieldfar-fieldclose)-dfdxf*g-2.*dfdxc*g);
d=(dfdxf*g+dfdxc*g-2*(fieldfar-fieldclose));
return fieldclose+dfdxc*p+c*pg*pg+d*pg*pg*pg;
}

//To minimize the numer of integrations, which have to be performed, this function estimates the four snapshots around the crossing of the world line and the observer's light cone.
int get_integration_range(struct light_rays_path *ray,double *rangemin,double *rangemax,double gridsize,int *SnapStart,int *SnapStop){
int i=-1,b,Start,Stop;
double dist,Mindist=-pow(gridsize*run.Ngrid,2.)*9.;
for (b=0;b<run.NSnap;++b) {// Find between which snapshots the world line and the light cone intersect
	dist=ray[b].dist*ray[b].dist-ray[b].TimeDistance*ray[b].TimeDistance;
	if ((dist>Mindist)&&(dist<0.0)) {Mindist=dist; i=b;}
}
if (i==-1) return 1;
if ((ray[i].dist-fabs(ray[i].vrHub) > run.maxdist+2*gridsize) || (ray[i].dist+fabs(ray[i].vrHub) < run.mindist-2*gridsize)) return 2;
Start=i-2;
Stop=i+2;
if (i>0) if (pow(ray[i-1].dist+ray[i-1].vrHub,2.)-pow(ray[i-1].TimeDistance-ray[i-1].vrHub,2.) < 0) --Start;
if (pow(ray[i].dist+ray[i].vrHub,2.)-pow(ray[i].TimeDistance-ray[i].vrHub,2.) > 0) ++Stop;
if (Start < 0) Start=0;
if (Stop > run.NSnap) Stop=run.NSnap;
SnapStart[0]=Start;
SnapStop[0]=Stop;
return 0;
}

//This function integrates the path, which was defined by "find_cells_along_ray" in light_rays.c
void integrate_path(struct light_rays_path *ray,double Hub,double *rangemin,double *rangemax,struct potential_container **Pot){
double MaxLineOfSight,IntegratedPot[12],IntStep[10],dr,Rcent,sum;
int i,j,k,n,l;

MaxLineOfSight=ray[0].PathIntersec[ray[0].CellsAlongPath-1];
ray[0].DRshift=ray[0].pot/(SPEEDLIGHT*Hub);
ray[0].DTshift=-ray[0].pot/(SPEEDLIGHT*Hub);
ray[0].Magni[0]=2.*ray[0].pot*SPEEDLIGHTSQUAREINV;
for (i=0;i<11;i++) IntegratedPot[i]=0.0;
for (i=1,n=0;(rangemin[i] < MaxLineOfSight)&&(i < run.NSnap);i++) {
	if (i==1) j=i-1; else j=i-2;
	if (i==run.NSnap-1) k=i; else k=i+1;
	for (;(n<ray[0].CellsAlongPath)&&(rangemax[i]>(ray[0].PathIntersec[n]+ray[0].PathIntersec[n-1])/2.0);n++) {
		l=ray[0].PathCell[n];
		if (n==0) dr=ray[0].PathIntersec[n]; else dr=ray[0].PathIntersec[n]-ray[0].PathIntersec[n-1];
		if (n==0) Rcent=ray[0].PathIntersec[n]*0.5; else Rcent=(ray[0].PathIntersec[n]+ray[0].PathIntersec[n-1])*0.5;

		IntegratedPot[0]+=simulation_time_interpol(Pot[j][l].val,Pot[i-1][l].val,Pot[i][l].val,Pot[k][l].val,rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;

		IntStep[0]=simulation_time_interpol(Pot[j][l].timeD,Pot[i-1][l].timeD,Pot[i][l].timeD,Pot[k][l].timeD,rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntegratedPot[1]+=IntStep[0];
		IntegratedPot[2]+=IntStep[0]*Rcent;

		IntStep[1]=simulation_time_interpol(Pot[j][l].firD[0],Pot[i-1][l].firD[0],Pot[i][l].firD[0],Pot[k][l].firD[0],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntegratedPot[3]+=IntStep[1];
		IntegratedPot[4]+=IntStep[1]*Rcent;

		//dist=simulation_time_interpol(dPhidy[j][l],dPhidy[i-1][l],dPhidy[i][l],dPhidy[k][l],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntStep[2]=simulation_time_interpol(Pot[j][l].firD[1],Pot[i-1][l].firD[1],Pot[i][l].firD[1],Pot[k][l].firD[1],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntegratedPot[5]+=IntStep[2];
		IntegratedPot[6]+=IntStep[2]*Rcent;

		IntStep[3]=simulation_time_interpol(Pot[j][l].firD[2],Pot[i-1][l].firD[2],Pot[i][l].firD[2],Pot[k][l].firD[2],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntegratedPot[7]+=IntStep[3];
		IntegratedPot[8]+=IntStep[3]*Rcent;

		IntStep[4]=simulation_time_interpol(Pot[j][l].secD[0],Pot[i-1][l].secD[0],Pot[i][l].secD[0],Pot[k][l].secD[0],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntStep[5]=simulation_time_interpol(Pot[j][l].secD[1],Pot[i-1][l].secD[1],Pot[i][l].secD[1],Pot[k][l].secD[1],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntStep[6]=simulation_time_interpol(Pot[j][l].secD[2],Pot[i-1][l].secD[2],Pot[i][l].secD[2],Pot[k][l].secD[2],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntStep[7]=simulation_time_interpol(Pot[j][l].croD[0],Pot[i-1][l].croD[0],Pot[i][l].croD[0],Pot[k][l].croD[0],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntStep[8]=simulation_time_interpol(Pot[j][l].croD[1],Pot[i-1][l].croD[1],Pot[i][l].croD[1],Pot[k][l].croD[1],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		IntStep[9]=simulation_time_interpol(Pot[j][l].croD[2],Pot[i-1][l].croD[2],Pot[i][l].croD[2],Pot[k][l].croD[2],rangemax[j],rangemax[i-1],rangemax[i],rangemax[k],Rcent)*dr;
		sum=IntStep[4]+IntStep[5]+IntStep[6];
		sum+=-2./MaxLineOfSight*(ray[0].normvec[0]*IntStep[1]+ray[0].normvec[1]*IntStep[2]+ray[0].normvec[2]*IntStep[3]);
		sum-=ray[0].normvec[0]*ray[0].normvec[0]*IntStep[4]+ray[0].normvec[1]*ray[0].normvec[1]*IntStep[5]+ray[0].normvec[2]*ray[0].normvec[2]*IntStep[6]+2.*ray[0].normvec[0]*ray[0].normvec[1]*IntStep[7]+2.*ray[0].normvec[0]*ray[0].normvec[2]*IntStep[8]+2.*ray[0].normvec[1]*ray[0].normvec[2]*IntStep[9];
		IntegratedPot[9]+=sum*Rcent;
		IntegratedPot[10]+=sum*Rcent*Rcent;
	}
}
//Now we compute the coordinate shift from the integrated quantities.
ray[0].DRshift+=4.0*IntegratedPot[0]*SPEEDLIGHTSQUAREINV;
ray[0].DRshift+=2.0*(1.0/(SPEEDLIGHT*Hub)+MaxLineOfSight*SPEEDLIGHTSQUAREINV)*IntegratedPot[1];
ray[0].DTshift-=2.0*IntegratedPot[1]/(SPEEDLIGHT*Hub);
ray[0].DRshift-=2.0*IntegratedPot[2]*SPEEDLIGHTSQUAREINV;

ray[0].D3shift[0]-=2.0*MaxLineOfSight*IntegratedPot[3]*SPEEDLIGHTSQUAREINV;
ray[0].D3shift[0]+=2.0*IntegratedPot[4]*SPEEDLIGHTSQUAREINV;

ray[0].D3shift[1]-=2.0*MaxLineOfSight*IntegratedPot[5]*SPEEDLIGHTSQUAREINV;
ray[0].D3shift[1]+=2.0*IntegratedPot[6]*SPEEDLIGHTSQUAREINV;

ray[0].D3shift[2]-=2.0*MaxLineOfSight*IntegratedPot[7]*SPEEDLIGHTSQUAREINV;
ray[0].D3shift[2]+=2.0*IntegratedPot[8]*SPEEDLIGHTSQUAREINV;

ray[0].Magni[0]+=-4.*IntegratedPot[0]/(SPEEDLIGHT*SPEEDLIGHT*MaxLineOfSight);
ray[0].Magni[0]+=-2.*(1.-SPEEDLIGHT/(Hub*MaxLineOfSight))*(-ray[0].pot*SPEEDLIGHTSQUAREINV+ray[0].vr/SPEEDLIGHT-2.*IntegratedPot[1]*SPEEDLIGHTSQUAREINV);
ray[0].Magni[0]+=2.*(IntegratedPot[9]-IntegratedPot[10]/MaxLineOfSight)*SPEEDLIGHTSQUAREINV;
ray[0].Magni[1]=2.*(IntegratedPot[9]-IntegratedPot[10]/MaxLineOfSight)*SPEEDLIGHTSQUAREINV;
ray[0].Magni[2]=-4.*IntegratedPot[0]/(SPEEDLIGHT*SPEEDLIGHT*MaxLineOfSight);
ray[0].Magni[3]=2.*ray[0].pot*SPEEDLIGHTSQUAREINV-2.*(1.-SPEEDLIGHT/(Hub*MaxLineOfSight))*(-ray[0].pot*SPEEDLIGHTSQUAREINV);
ray[0].Magni[4]=-2.*(1.-SPEEDLIGHT/(Hub*MaxLineOfSight))*(ray[0].vr/SPEEDLIGHT);
ray[0].Magni[5]=-2.*(1.-SPEEDLIGHT/(Hub*MaxLineOfSight))*(-2.*IntegratedPot[1]*SPEEDLIGHTSQUAREINV);

}

