.PHONY: all clean clean-all lulesh libcm redis flann

all: lulesh

lulesh: libcm
	${MAKE} -C LULESH

libcm: redis flann
	${MAKE} -C CM/exec REDIS=yes FLANN=yes

redis:
	${MAKE} -C redis

flann:
	${MAKE} -C flann

clean:
	${MAKE} -C CM/exec clean
	${MAKE} -C LULESH clean

clean-all: clean
	${MAKE} -C redis clean
	${MAKE} -C flann clean
