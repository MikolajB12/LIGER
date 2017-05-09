EXEC   = Liger

OBJS   = main.o src/ParaIO.o src/integration.o src/light_rays.o src/load_snap.o src/load_snap_gdt.o src/potential.o

INCL   = src/vardef.h Makefile


OPT   += -DMEASURETIME

#PAR   += -DNOOPENMP


OPTIONS =  $(OPT) $(PAR)

#FFTW_INCL = -I/users/aludlow/InstallPrograms/fftw-2.1.5/include/
#FFTW_LIBS = -L/users/aludlow/InstallPrograms/fftw-2.1.5/lib/

CC       =  gcc  #mpicc

OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)


FFTW_LIB =  -lfftw3_threads  -lfftw3
#ifeq (SINGLEFFTW,$(findstring SINGLEFFTW,$(OPT)))
#  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
#else
#  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
#endif

LIBS   =   -lm  $(FFTW_LIB)

CFLAGS =   $(OPTIONS)  $(OPTIMIZE)

ifneq ($(PAR),-DNOOPENMP)
  CFLAGS += -fopenmp
  LIBS += -fopenmp
endif

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC)



