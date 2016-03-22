# flags for gcc
#ifeq ($(CC), gcc)
FC = gfortran 
CXX = g++
FFLAGS = -ffree-line-length-none -O3 -fdefault-real-8
FLIBS = -lm -lgfortran -lquadmath
#endif
# flags for intel
ifeq ($(COMPILER), intel)
FC = ifort
CXX = icc
FFLAGS = -real-size 64
FLIBS = -lm -lifcore
endif


FILES = vpEtxt_interface_base.f90 vpscVpETxt_interface.f90 libraryR8.f vpsc_as_N2O.f init.f90 TaylorFort.f90

SRCS = $(wildcard $(FILES))

OBJS = vpEtxt_interface_base.o vpscVpETxt_interface.o libraryR8.o vpsc_as_N2O.o init.o
MODS = vpetxt_interface_base.mod vpscvpetxt_interface.mod

Taylor: TaylorFort.f90
	$(FC) -c TaylorFort.f90

modules: $(MODS) $(OBJS)

vpetxt_interface_base.mod vpEtxt_interface_base.o: vpEtxt_interface_base.f90
	$(FC) $(FFLAGS) -c $<

vpscvpetxt_interface.mod vpscVpETxt_interface.o: vpscVpETxt_interface.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.mod *.o

print-%  : ; @echo $* = $($*)
