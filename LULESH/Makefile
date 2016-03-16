###### CONFIGURATION START ############

# Use MPI
COEVP_MPI =

# Silo (www.hdfgroup.org/HDF5) and HDF5 (www.hdfgroup.org/HDF5) library
# locations.  These macros can be left undefined if visualization is not required.
SILO_LOC = 

# HDF5 library location (www.hdfgroup.org/HDF5).  This macro can
# be undefined if visualization is not required.
HDF5_LOC = 

# BLAS and LAPACK library locations.  It may not be necessary
# to specify these depending upon the linker defaults
LAPACK_LOC =
BLAS_LOC =

# FLANN library location
FLANN_LOC =

# Redis library location
REDIS_LOC =

LOGGER_LOC=

#CXX = g++
#CXXFLAGS = -std=c++0x -g -O3 -fopenmp
CXXFLAGS = -std=c++11 -g -O3

###### CONFIGURATION END ############

CM_LIB = ../CM/lib/libcm.a

INCLUDE_DIRS = -I. -I../CM/include

ifeq ($(COEVP_MPI), yes)
CXXFLAGS += -DCOEVP_MPI
CXX= mpic++
endif

ifeq ($(CHARM), yes)
CXXFLAGS += -D__CHARMC__
ifeq ($(COEVP_MPI), yes)
$(error You cannot combine charm with mpi)
endif
endif

ifneq ($(strip $(FLANN_LOC)),)
INCLUDE_DIRS +=  -I$(FLANN_LOC) -I$(FLANN_LOC)/flann
CXXFLAGS += -DFLANN
endif

ifneq ($(strip $(REDIS_LOC)),)
INCLUDE_DIRS +=  -I$(REDIS_LOC)
CXXFLAGS += -DREDIS
endif

SILO_LIB =
ifneq ($(strip $(SILO_LOC)),)
INCLUDE_DIRS += -I$(SILO_LOC)/include/
SILO_LIB += -Wl,-rpath,$(realpath $(SILO_LOC))/lib -L $(SILO_LOC)/lib -lsiloh5
CXXFLAGS += -DSILO 
endif


LOG_LIB =
LOGGER_LIB =
ifneq ($(strip $(LOGGER_LOC)),)
INCLUDE_DIRS += -I$(LOGGER_LOC)
LOGGER_LIB += -Wl,-rpath,$(realpath $(LOGGER_LOC)) -L $(LOGGER_LOC) -llogger
CXXFLAGS += -DLOGGER
LOG_LIB = ../logger/liblogger.a     # make lulesh depend on logger library
endif

PROTO_LIB =
PROTOBUF_LIB =
PBSRCS   = $(wildcard *.proto)
PBOBJS   = $(PBSRCS:.proto=.pb.o)
PBGENS   = $(PBSRCS:.proto=.pb.cc) $(PBSRCS:.proto=.pb.h)
ifneq ($(strip $(PROTOBUF_LOC)),)
INCLUDE_DIRS += -I$(PROTOBUF_LOC)/include
PROTOBUF_LIB += -Wl,-rpath,$(realpath $(PROTOBUF_LOC))/lib -L $(PROTOBUF_LOC)/lib -lprotobuf
CXXFLAGS += -DPROTOBUF
PROTO_LIB = $(PROTOBUF_LOC)/lib/libprotobuf.a   # make lulesh depend on protobuf library
PROTOC = $(PROTOBUF_LOC)/bin/protoc
endif

HDF5_LIB =
ifneq ($(strip $(HDF5_LOC)),)
HDF5_LIB += -Wl,-rpath,$(HDF5_LOC)/lib -L $(HDF5_LOC)/lib -lhdf5
endif

LAPACK = $(shell pkg-config --silence-errors --libs lapack blas || true)
ifeq ($(strip $(LAPACK)),)
ifneq ($(strip $(LAPACK_LOC)),)
LAPACK += -L$(LAPACK_LOC)
endif
LAPACK += -llapack

ifneq ($(strip $(BLAS_LOC)),)
LAPACK += -L$(BLAS_LOC)
endif
LAPACK += -lblas
endif

ifeq ($(FSTRACE), yes)
CXXFLAGS += -DFSTRACE
endif

CXXFLAGS += $(INCLUDE_DIRS)

OPTFLAGS = $(CXXFLAGS) $(CCFLAGS)

all: lulesh

dummy: ;

.cxxflags: dummy
	@[ -f $@ ] || touch $@
	@echo "CXXFLAGS=$(CXXFLAGS)" | cmp -s $@ - || echo "CXXFLAGS=$(CXXFLAGS)" > $@

# The Makefile doesn't work correctly.
# If the .pb.cc and .pb.h file DO NOT exist, they will NOT be built by this target.
# If they do exist, and the .proto file is newer, the target WILL be run.
# To get the makefile to work properly, need to find out how to run protoc when the
# .pb.cc and .pb.h file DO NOT PREVIOUSLY EXIST.
# So, for now, before forst execution of Makefile, run (from within LULESH directory):
# ../serverize/protobuf/bin/protoc --cpp_out=. advance.proto
%.pb.cc %.pb.h: %.proto
	$(PROTOC) --cpp_out=. $<

%.o: %.cc .cxxflags
	$(CXX) -c $(CXXFLAGS) -c -o $@ $< 

LULESH_SRC=$(filter-out $(PBGENS) lctest.cc,$(shell ls *.cc))
LULESH_OBJS=$(LULESH_SRC:.cc=.o)
ifneq ($(strip $(PROTOBUF_LOC)),)
LULESH_OBJS+=$(PBOBJS)
shims.o: advance.pb.h; 
endif

# Al added the -lrt library on 2/22/2016 because Travis builds were failing with an undefined reference
# to clock_gettime. Documentation says that -lrt is required, but it was working OK without on Darwin.
# It does build successfully on Darwin with -lrt too.
lulesh:  $(LULESH_OBJS) $(CM_LIB) $(LOG_LIB)
	$(CXX) $(OPTFLAGS) -o $@ $(LULESH_OBJS) $(CM_LIB) $(SILO_LIB) $(HDF5_LIB) $(LOGGER_LIB) $(PROTOBUF_LIB) $(LAPACK) -lrt $(EXTRA_LIBS) -lm

clean:
	rm -f *~ lulesh $(LULESH_OBJS) $(PBGENS) $(PBOBJS)
