#
# ===== Build options begin here =====
#

#
# Specify the C++ compiler.  This is also used for linking the executable.
#

CXX = g++

#
# Specify compiler options (optimization level, debugging options, etc.)
#

CXXFLAGS = -O3

#
# Specify the location of the BLAS and LAPACK libraries.  Depending upon the
# directories searched by the linker by default, it may not be necessary to
# specify these, in which case they can be left undefined (empty).  You can first
# try building with these macros undefined (empty) to find out.  If the libraries
# aren't found, or if a specific implementation of BLAS/LAPACK located elsewhere
# is desired, set these macros to the location of libblas.a and liblapack.a,
# respectively.
#

LAPACK_LOC =
BLAS_LOC =

#
# The M-tree database can be read/written from/to a file in HDF5 format, but
# this is optional.  If such functionality is desired, then specify the location
# of the HDF5 library (i.e., the top level directory containing the bin, include
# and lib subdirectories).  Otherwise just leave this macro undefined (empty).
#

HDF5_LOC =
HDF5_DEFINES = -DH5_USE_16_API

#
# ===== Build options end here =====
#

LIBS =
ifneq ($(strip $(LAPACK_LOC)),)
LIBS += -L$(LAPACK_LOC)
endif
LIBS += -llapack

ifneq ($(strip $(BLAS_LOC)),)
LIBS += -L$(BLAS_LOC)
endif
LIBS += -lblas

ifeq ($(strip $(HDF5_LOC)),)
CXXFLAGS += -UHAVE_PKG_hdf5
else
CXXFLAGS += -DHAVE_PKG_hdf5 -I$(HDF5_LOC)/include $(HDF5_DEFINES)
LIBS += -L$(HDF5_LOC)/lib -lhdf5
endif

# Tensors

TENSOR_LOC  = ../src/tensors
TENSOR_SRC  = $(TENSOR_LOC)/tensor.cc $(TENSOR_LOC)/xtensor.cc $(TENSOR_LOC)/tutils.cc
TENSOR_INCL = $(TENSOR_SRC:.cc=.h)
TENSOR_OBJS = $(TENSOR_SRC:.cc=.o)
TENSOR_DEPS = $(TENSOR_SRC:.cc=.d)

# Adaptive sampling

AS_LOC       = ../src/adaptive_sampling
INTERPDB_LOC = $(AS_LOC)/interpolation_database
INTERPDB_SRC = $(shell ls $(INTERPDB_LOC)/base/*.cc $(INTERPDB_LOC)/key_db/*.cc)
INTERP_LOC   = $(AS_LOC)/interpolation
INTERP_SRC   = $(shell ls $(INTERP_LOC)/*/*.cc)
DB_LOC       = $(AS_LOC)/database
DB_SRC       = $(shell ls $(DB_LOC)/base/*.cc $(DB_LOC)/mtreedb/*.cc)
UTILS_LOC    = $(AS_LOC)/utils
UTILS_SRC    = $(shell ls $(UTILS_LOC)/toolbox/*/*.cc)
MTL_LOC      = $(AS_LOC)/mtl_headers
AS_SRC       = $(AS_LOC)/AdaptiveSampler.cc $(INTERPDB_SRC) $(INTERP_SRC) $(DB_SRC) $(UTILS_SRC)
AS_OBJS      = $(AS_SRC:.cc=.o)
AS_DEPS      = $(AS_SRC:.cc=.d)

# Constitutive models

CM_LOC   = ../src/constitutive_models
CM_SRC   = $(shell ls $(CM_LOC)/*.cc)
CM_INCL  = $(CM_SRC:.cc=.h)
CM_OBJS  = $(CM_SRC:.cc=.o)
CM_DEPS  = $(CM_SRC:.cc=.d)

# EOS models

EOS_LOC  = ../src/eos_models
EOS_SRC  = $(shell ls $(EOS_LOC)/*.cc)
EOS_INCL = $(EOS_SRC:.cc=.h) $(EOS_LOC)/EOS.h
EOS_OBJS = $(EOS_SRC:.cc=.o)
EOS_DEPS = $(EOS_SRC:.cc=.d)

# Fine scale models

FS_LOC  = ../src/fine_scale_models
FS_SRC  = $(shell ls $(FS_LOC)/*.cc)
FS_INCL = $(FS_SRC:.cc=.h) $(FS_LOC)/FineScale.h
FS_OBJS = $(FS_SRC:.cc=.o)
FS_DEPS = $(FS_SRC:.cc=.d)

# Standalone kinematic system test

KINTEST_OBJS = kintest.o $(CM_OBJS) $(EOS_OBJS) $(FS_OBJS) $(AS_OBJS) $(TENSOR_OBJS)
KINTEST_DEPS = kintest.d

KINTEST2_OBJS = kintest2.o $(CM_OBJS) $(EOS_OBJS) $(FS_OBJS) $(AS_OBJS) $(TENSOR_OBJS)
KINTEST2_DEPS = kintest2.d

# Constitutive model library

CMLIB_SRC  = $(CM_SRC) $(EOS_SRC) $(FS_SRC) $(AS_SRC) $(TENSOR_SRC)
CMLIB_OBJS = $(CM_OBJS) $(EOS_OBJS) $(FS_OBJS) $(AS_OBJS) $(TENSOR_OBJS)
CMLIB_INCL = $(CM_INCL) $(EOS_INCL) $(FS_INCL) $(TENSOR_INCL)

CXXFLAGS += -std=c++0x
CXXFLAGS += -I$(AS_LOC) -I$(MTL_LOC) -I$(INTERPDB_LOC) -I$(INTERP_LOC) -I$(DB_LOC) -I$(UTILS_LOC)
CXXFLAGS += -I$(TENSOR_LOC) -I$(CM_LOC) -I$(EOS_LOC) -I$(FS_LOC)
CXXFLAGS += -DDBL_SNAN_IS_BROKEN -DFLT_SNAN_IS_BROKEN

MAKEDEPEND = gcc -MM -MT $*.o $(CXXFLAGS) -o $*.d.tmp $<

default: libcm

kintest: $(KINTEST_OBJS)
	$(CXX) $(CXXFLAGS) $(KINTEST_OBJS) $(LIBS) -o kintest

kintest2: $(KINTEST2_OBJS)
	$(CXX) $(CXXFLAGS) $(KINTEST2_OBJS) $(LIBS) -o kintest2

libcm:  $(CMLIB_OBJS)
	ar crs libcm.a $(CMLIB_OBJS)
	cp $(CMLIB_INCL) ../include
	mv libcm.a ../lib

%.o : %.cc
	@$(MAKEDEPEND)
	cp $*.d.tmp $*.d
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' -e '/^$$/ d' -e 's/$$/ :/' < $*.d.tmp >> $*.d
	rm -f $*.d.tmp
	$(CXX) -c $(CXXFLAGS) -c -o $@ $< 

-include $(CM_DEPS)
-include $(AS_DEPS)
-include $(FS_DEPS)
-include $(TENSOR_DEPS)
-include $(EOS_DEPS)
-include $(KINTEST_DEPS)

clean:
	$(RM) $(AS_OBJS) $(AS_DEPS) $(CM_OBJS) $(CM_DEPS) $(FS_OBJS) $(FS_DEPS)
	$(RM) $(TENSOR_OBJS) $(TENSOR_DEPS) $(EOS_OBJS) $(EOS_DEPS) $(KINTEST_OBJS) $(KINTEST_DEPS) *.a

realclean: clean
	$(RM) kintest ../include/* ../lib/*


