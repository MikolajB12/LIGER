/*  This file is part of Liger, see main.c and the User guide.
    Copyright (C) 2017  Mikolaj Borzyszkowski, Argelander Institut für Astronomie, University of Bonn */

#ifndef _POTENTIAL_H
#define _POTENTIAL_H
#include "vardef.h"
void fill_potential_container(float *PotInput,int Ngrid,double boxsize,struct potential_container *PotCont); //Computes all spacial derivatives.
int overdens(double *Alldens,int Ngrid);
double *density_field(struct particle_data *P,struct header_snapshot header1,unsigned long Ngrid,int flag);
double weight(struct particle_data part,double gridx,double gridy,double gridz,double gridsize,int flag);
double get_potential_value(double posx,double posy,double posz,double gridsize,int Ngrid,struct potential_container *Pot);
float *grav_potential(double **dens,int Ngrid,struct header_snapshot header1);
void compute_potential_timeder(int Ngrid3,int Nsnap,double *rangemax,struct header_snapshot *Head,struct potential_container **PotCont);
int shift_koor(int a,int b,int c,int Ngrid);
double get_distance_square_periodic(double x1,double x2,double x3,double y1,double y2,double y3,double boxsize);
double PeriDist(double dist,double box);
void get_potential(char *fname,float **pot,int Ngrid,struct header_snapshot Head,char *path);
#endif
