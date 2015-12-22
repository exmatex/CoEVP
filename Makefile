.PHONY: all clean clean-all lulesh libcm redis flann twemproxy silo test

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
TWEMPROXY=yes
ifeq ($(TWEMPROXY),yes)
libcm: twemproxy
endif

lulesh: LULESH/lulesh

LULESH/lulesh: libcm
	${MAKE} -C LULESH FLANN_LOC=$(FLANN_LOC) SILO_LOC=$(SILO_LOC) REDIS_LOC=$(REDIS_LOC)

libcm:
	${MAKE} -C CM/exec REDIS=$(REDIS) FLANN=$(FLANN) TWEMPROXY=$(TWEMPROXY)

redis:
	${MAKE} -C redis

silo:
	${MAKE} -C silo

flann:
	${MAKE} -C flann

twemproxy:
	${MAKE} -C twemproxy

clean:
	${MAKE} -C CM/exec clean
	${MAKE} -C LULESH clean
	rm -rf test/*.silo

clean-all: clean
	${MAKE} -C redis clean
	${MAKE} -C flann clean
	${MAKE} -C silo clean
	${MAKE} -C twemproxy clean

reference: LULESH/lulesh
	@[ "$(SILO)" = "yes" ] || { echo "make test needs SILO=yes" && exit 1; }
	mkdir -p test/reference
	cd test/reference && ../../LULESH/lulesh

STEPS=500
#bit hackish, but let's assume we have $(STEPS) steps
test/taylor_$(STEPS).silo: LULESH/lulesh
	@[ "$(SILO)" = "yes" ] || { echo "make test needs SILO=yes" && exit 1; }
	mkdir -p test
	cd test && ../LULESH/lulesh

test: test/taylor_$(STEPS).silo silo
	@[ -x "$(SILODIFF)" ] || { echo "SILODIFF=$(SILODIFF) seems to be wrong" && exit 1; }
	$(SILODIFF) test/reference test > test/diff
	@[ ! -s test/diff ] || { echo "Difference in files" && exit 1; }
