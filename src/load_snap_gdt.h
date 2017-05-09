/*  This file is part of Liger, see main.c and the User guide.
    Copyright (C) 2017  Mikolaj Borzyszkowski, Argelander Institut für Astronomie, University of Bonn */

#ifndef _LOAD_SNAP_GDT_H
#define _LOAD_SNAP_GDT_H
#include "vardef.h"
struct header_snapshot load_header_gdt(char *fname);
int load_snapshot_gdt(char *fname, int files,struct particle_data **Pptr);
int load_snapshot_gdt_partial(char *fname,int files,int NStart,int NSelect,struct particle_data **Pptr);
int allocate_memory(int NumPart,struct particle_data **P);
int allocate_memory_partial(int NSelection,struct particle_data **P);
#endif
