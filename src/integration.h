/*  This file is part of Liger, see main.c and the User guide.
    Copyright (C) 2017  Mikolaj Borzyszkowski, Argelander Institut für Astronomie, University of Bonn */

#ifndef _INTEGRATION_H
#define _INTEGRATION_H
#include "vardef.h"
double Hub(double a,double Om,double OL);
double calc_comoving_distance(double red,double dz,double Om,double OL);
double simulation_time_interpol(double fieldcloser,double fieldclose,double fieldfar,double fieldfarer,double distcloser,double distclose,double distfar,double distfarer,double dist);
int get_integration_range(struct light_rays_path *ray,double *rangemin,double *rangemax,double gridsize,int *SnapStart,int *SnapStop);
void integrate_path(struct light_rays_path *ray,double Hub,double *rangemin,double *rangemax,struct potential_container **Pot);
#endif
