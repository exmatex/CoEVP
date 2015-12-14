.PHONY: all clean clean-all lulesh libcm redis flann silo

all: lulesh

FLANN=yes
ifeq ($(FLANN),yes)
FLANN_LOC=../flann/flann/src/cpp
libcm: flann
endif
REDIS=yes
ifeq ($(REDIS),yes)
libcm: redis
endif
SILO=yes
ifeq ($(SILO),yes)
SILO_LOC=../silo/silo
libcm: silo
endif

lulesh: libcm
	${MAKE} -C LULESH FLANN_LOC=$(FLANN_LOC) SILO_LOC=$(SILO_LOC)

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

clean-all: clean
	${MAKE} -C redis clean
	${MAKE} -C flann clean
	${MAKE} -C silo clean
