# flags for gcc
#ifeq ($(CC), gcc)
FC ?= gfortran 
CXX ?= g++
#FFLAGS = -ffree-line-length-none -O3 -fdefault-real-8
FFLAGS = -ffree-line-length-none -O3 -g -fdefault-real-8
FLIBS = -lm -lgfortran -lquadmath
#endif
# flags for intel
ifeq ($(COMPILER), intel)
FC = ifort
CXX = icc
FFLAGS = -real-size 64
FLIBS = -lm -lifcore
endif


FILES = vpetxt_interface_base.f90 vpscvpetxt_interface.f90 libraryR8.f vpsc_as_N2O.f init.f90 TaylorFort.f90

SRCS = $(wildcard $(FILES))

OBJS = vpetxt_interface_base.o vpscvpetxt_interface.o libraryR8.o vpsc_as_N2O.o init.o
MODS = vpetxt_interface_base.mod vpscvpetxt_interface.mod

all: Taylor modules

Taylor: TaylorFort.o

modules: $(MODS) $(OBJS)

vpscvpetxt_interface.o: vpetxt_interface_base.mod
init.o: vpscvpetxt_interface.mod
caller.o: vpscvpetxt_interface.mod

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

%.mod: %.o;

clean:
	rm -f *.mod *.o

print-%  : ; @echo $* = $($*)
