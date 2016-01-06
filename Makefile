.PHONY: all clean clean-all lulesh libcm redis flann silo test

all: lulesh

FLANN=yes
ifeq ($(FLANN),yes)
FLANN_LOC=../flann/flann/src/cpp
libcm: flann
endif
REDIS=yes
ifeq ($(REDIS),yes)
REDIS_LOC=../redis/hiredis
libcm: redis
endif
SILO=yes
ifeq ($(SILO),yes)
SILO_LOC=../silo/silo
SILODIFF=silo/silo/bin/silodiff
libcm: silo
endif

lulesh: LULESH/lulesh

LULESH/lulesh: libcm
	${MAKE} -C LULESH FLANN_LOC=$(FLANN_LOC) SILO_LOC=$(SILO_LOC) REDIS_LOC=$(REDIS_LOC)

libcm:
	${MAKE} -C CM/exec REDIS=$(REDIS) FLANN=$(FLANN)

redis:
	${MAKE} -C redis

silo:
	${MAKE} -C silo

flann:
	${MAKE} -C flann

clean:
	${MAKE} -C CM/exec clean
	${MAKE} -C LULESH clean
	rm -rf test/*.silo

clean-all: clean
	${MAKE} -C redis clean
	${MAKE} -C flann clean
	${MAKE} -C silo clean

get_reference:
	mkdir -p test/reference
	git clone https://github.com/exmatex/CoEVP_reference.git test/reference

LULESH_OPTS=-p 4
reference: LULESH/lulesh
	@[ "$(SILO)" = "yes" ] || { echo "make test needs SILO=yes" && exit 1; }
	mkdir -p test/reference
	cd test/reference && ../../LULESH/lulesh $(LULESH_OPTS)
	rm test/reference/*.silo #remove total files

dummy: ;

test/.mpirunflags: dummy
	@[ -f $@ ] || touch $@
	@echo "MPIRUN=$(MPIRUN)" | cmp -s $@ - || echo "MPIRUN=$(MPIRUN)" > $@

STEPS=0500
#bit hackish, but let's assume we have $(STEPS) steps
test/taylor_$(STEPS).silo: LULESH/lulesh test/.mpirunflags
	@[ "$(SILO)" = "yes" ] || { echo "make test needs SILO=yes" && exit 1; }
	mkdir -p test
	cd test && $(MPIRUN) ../LULESH/lulesh $(LULESH_OPTS)

SILODIFF_OPTS=-A 1e-8 -E _hdf5libinfo
test: test/taylor_$(STEPS).silo silo
	@[ -x "$(SILODIFF)" ] || { echo "SILODIFF=$(SILODIFF) seems to be wrong" && exit 1; }
	$(SILODIFF) ${SILODIFF_OPTS} test/reference test > test/diff
	@[ ! -s test/diff ] || { echo "Difference in files" && head -n 50 test/diff && exit 1; }
