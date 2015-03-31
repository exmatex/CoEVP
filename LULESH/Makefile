###### CONFIGURATION START ############

# Silo (www.hdfgroup.org/HDF5) and HDF5 (www.hdfgroup.org/HDF5) library
# locations.  These macros can be left undefined if visualization is not required.

#SILO_LOC = /usr/local/tools/silo-4.8
SILO_LOC = 

# HDF5 library location (www.hdfgroup.org/HDF5).  This macro can
# be undefined if visualization is not required.
#HDF5_LOC = /usr/local/tools/hdf5-intel-serial-1.8.8
HDF5_LOC = 

# BLAS and LAPACK library locations.  It may not be necessary
# to specify these depending upon the linker defaults
LAPACK_LOC =
BLAS_LOC =

CXX = g++
CXXFLAGS = -g -O3 -fopenmp

###### CONFIGURATION END ############

CM_LIB = ../CM/lib/libcm.a

INCLUDE_DIRS = -I../CM/include

SILO_LIB =
ifneq ($(strip $(SILO_LOC)),)
INCLUDE_DIRS += -I$(SILO_LOC)/include/
SILO_LIB += -L $(SILO_LOC)/lib -lsiloh5 
endif

HDF5_LIB =
ifneq ($(strip $(HDF5_LOC)),)
HDF5_LIB += -Wl,-rpath,$(HDF5_LOC)/lib -L $(HDF5_LOC)/lib -lhdf5
endif

LAPACK =
ifneq ($(strip $(LAPACK_LOC)),)
LAPACK += -L$(LAPACK_LOC)
endif
LAPACK += -llapack

ifneq ($(strip $(BLAS_LOC)),)
LAPACK += -L$(BLAS_LOC)
endif
LAPACK += -lblas

lulesh: lulesh.cc $(CM_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) $< $(CM_LIB) $(SILO_LIB) $(HDF5_LIB) $(LAPACK) -lm -o $@

clean:
	rm -f *~ lulesh
