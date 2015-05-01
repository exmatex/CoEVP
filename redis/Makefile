REDIS_PV=2.8.19
HIREDIS_PV=0.13.0

REDIS_SERVER=redis-${REDIS_PV}/src/redis-server
HIREDIS_LIB=hiredis-${HIREDIS_PV}/libhiredis.a

all: ${REDIS_SERVER} $(HIREDIS_LIB)

redis-${REDIS_PV}.tar.gz:
	wget -O $@ https://github.com/antirez/redis/archive/${REDIS_PV}.tar.gz

hiredis-${HIREDIS_PV}.tar.gz:
	wget -O $@ https://github.com/redis/hiredis/archive/v${HIREDIS_PV}.tar.gz

$(REDIS_SERVER): redis-${REDIS_PV}.tar.gz
	tar -xzvf $<
	make -C redis-${REDIS_PV}
	echo Please use $(REDIS_SERVER)

$(HIREDIS_LIB): hiredis-${HIREDIS_PV}.tar.gz
	tar -xzvf $<
	make -C hiredis-${HIREDIS_PV}

clean:
	rm -rf redis-${REDIS_PV}
	rm -rf hiredis-${HIREDIS_PV}

clean-tar:
		rm -rf *.tar.gz
