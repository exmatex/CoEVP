!FFLAGS = -gdwarf-2 -ffree-line-length-0 -pg
!FFLAGS = -gdwarf-2 -pg -check all
!FFLAGS = -gdwarf-2 -check all -real-size 64
# flags for intel
#ifeq ($(CC), icc)
FC = ifort
CXX = icc
FFLAGS = -real-size 64
FLIBS = -lm -lifcore
#endif
# flags for gcc
#ifeq ($(CC), gcc)
#FC = gfortran 
#CXX = g++
#FFLAGS = -ffree-line-length-none -O3 -fdefault-real-8
#FLIBS = -lm -lgfortran -lquadmath
#endif


FILES = vpEtxt_interface_base.f90 vpscVpETxt_interface.f90 libraryR8.f vpsc_as_N2O.f init.f90 TaylorFort.f90

SRCS = $(wildcard $(FILES))

OBJS = vpEtxt_interface_base.o vpscVpETxt_interface.o libraryR8.o vpsc_as_N2O.o init.o
MODS = vpETxt_interface_base.mod vpscVpETxt_interface.mod

Taylor: TaylorFort.f90
	$(FC) -c TaylorFort.f90

modules: $(MODS) $(OBJS) $(SRCS)
	$(FC) $(FFLAGS) -c $(SRCS)

vpETxt_interface_base.mod: vpEtxt_interface_base.f90
	$(FC) $(FFLAGS) -c $<

vpscVpETxt_interface.mod: vpscVpETxt_interface.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm *.mod *.o

print-%  : ; @echo $* = $($*)