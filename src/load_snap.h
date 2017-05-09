/*  This file is part of Liger, see main.c and the User guide.
    Copyright (C) 2017  Mikolaj Borzyszkowski, Argelander Institut für Astronomie, University of Bonn */

#ifndef _LOAD_SNAP_H
#define _LOAD_SNAP_H
#include "vardef.h"
struct header_snapshot load_snapshot_header(int FileNr);
struct header_snapshot load_galaxies_header(int FileNr);
int load_snapshot(int FileNr, int files,struct particle_data **Pptr);
int load_galaxies_partial(int FileNr,int files,int NStart,int NSelect,struct particle_data **Pptr);//Currently for Gadget only files=1 is possible.
#endif
