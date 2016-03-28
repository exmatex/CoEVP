.PHONY: all clean clean-all lulesh libcm redis flann silo test logger twemproxy protobuf circle

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
# LOGGER depends on REDIS=yes
# IF REDIS=no, logging ends up being a NO OP
LOGGER=no
ifeq ($(LOGGER)$(REDIS), yesno)
LOGGER=no
endif
ifeq ($(LOGGER)$(REDIS), yesyes)
LOGGER_LOC=../logger
libcm: logger
endif
PROTOBUF=no
ifeq ($(PROTOBUF),yes)
PROTOBUF_LOC=../serverize/protobuf
CIRCLE_LOC=../serverize/circle
libcm: protobuf
endif
FSTRACE=no
USE_SSL=yes
ifeq ($(USE_SSL),no)
CURLFLAG=-k
WGETFLAG=--no-check-certificate
endif

lulesh: LULESH/lulesh

LULESH/lulesh: libcm
	${MAKE} -C LULESH FLANN_LOC=$(FLANN_LOC) SILO_LOC=$(SILO_LOC) REDIS_LOC=$(REDIS_LOC) LOGGER_LOC=$(LOGGER_LOC) PROTOBUF_LOC=$(PROTOBUF_LOC) CIRCLE_LOC=$(CIRCLE_LOC) FSTRACE=$(FSTRACE) 

libcm:
	${MAKE} -C CM/exec REDIS=$(REDIS) FLANN=$(FLANN) TWEMPROXY=$(TWEMPROXY) FSTRACE=$(FSTRACE) LOGGER=$(LOGGER) PROTOBUF=$(PROTOBUF)

redis:
	${MAKE} -C redis CURLFLAG=$(CURLFLAG) WGETFLAG=$(WGETFLAG)

silo:
	${MAKE} -C silo CURLFLAG=$(CURLFLAG) WGETFLAG=$(WGETFLAG)

flann:
	${MAKE} -C flann CURLFLAG=$(CURLFLAG) WGETFLAG=$(WGETFLAG)

twemproxy:
	${MAKE} -C twemproxy CURLFLAG=$(CURLFLAG) WGETFLAG=$(WGETFLAG)

logger: redis
	${MAKE} -C logger REDIS=$(REDIS) REDIS_LOC=$(REDIS_LOC)

protobuf:
	${MAKE} -C serverize protobuf CURLFLAG=$(CURLFLAG) WGETFLAG=$(WGETFLAG)

circle:
	${MAKE} -C serverize/circle

clean:
	${MAKE} -C CM/exec realclean
	${MAKE} -C LULESH clean
	rm -rf test/*.silo

clean-all: clean
	${MAKE} -C redis clean
	${MAKE} -C flann clean
	${MAKE} -C silo clean
	${MAKE} -C twemproxy clean
	${MAKE} -C logger clean
	${MAKE} -C serverize/protobuf clean
	${MAKE} -C serverize/circle clean

get_reference:
	mkdir -p test/reference
	git clone https://github.com/exmatex/CoEVP_reference.git test/reference

LULESH_OPTS=-p 4 -v 20
reference: LULESH/lulesh
	@[ "$(SILO)" = "yes" ] || { echo "make test needs SILO=yes" && exit 1; }
	mkdir -p test/reference
	cd test/reference && ../../LULESH/lulesh $(LULESH_OPTS)
	rm test/reference/*.silo #remove total files

dummy: ;

test/.mpirunflags: dummy
	mkdir -p test
	@[ -f $@ ] || touch $@
	@echo "MPIRUN=$(MPIRUN)" | cmp -s $@ - || echo "MPIRUN=$(MPIRUN)" > $@

test/.luleshopts: dummy
	mkdir -p test
	@[ -f $@ ] || touch $@
	@echo "LULESH_OPTS=$(LULESH_OPTS)" | cmp -s $@ - || echo "LULESH_OPTS=$(LULESH_OPTS)" > $@

STEPS=0500
#bit hackish, but let's assume we have $(STEPS) steps
test/taylor_$(STEPS).silo: LULESH/lulesh test/.mpirunflags test/.luleshopts
	mkdir -p test
	ifdef $(NTASKS)
		cp LULESH/taskhandler test/.
	endif
	cd test && NTASKS=$(NTASKS) NHANDLERS=$(NHANDLERS) $(MPIRUN) ../LULESH/lulesh $(LULESH_OPTS)

SILODIFF_OPTS=-A 1e-8 -E _hdf5libinfo
test: test/taylor_$(STEPS).silo silo
	@[ -x "$(SILODIFF)" ] || { echo "SILODIFF=$(SILODIFF) seems to be wrong" && exit 1; }
	$(SILODIFF) ${SILODIFF_OPTS} test/reference test > test/diff
	@[ ! -s test/diff ] || { echo "Difference in files" && head -n 50 test/diff && exit 1; }
