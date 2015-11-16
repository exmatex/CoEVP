REDIS_PV=2.8.22
HIREDIS_PV=0.13.3

REDIS_SERVER=redis-${REDIS_PV}/src/redis-server
HIREDIS_LIB=hiredis-${HIREDIS_PV}/libhiredis.a

all: ${REDIS_SERVER} $(HIREDIS_LIB)

redis-${REDIS_PV}.tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://github.com/antirez/redis/archive/${REDIS_PV}.tar.gz
else
	curl -o $@ -L https://github.com/antirez/redis/archive/${REDIS_PV}.tar.gz
endif

hiredis-${HIREDIS_PV}.tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://github.com/redis/hiredis/archive/v${HIREDIS_PV}.tar.gz
else
	curl -o $@ -L https://github.com/redis/hiredis/archive/v${HIREDIS_PV}.tar.gz
endif

$(REDIS_SERVER): redis-${REDIS_PV}.tar.gz
	tar -xzvf $<
	make -C redis-${REDIS_PV}
	ln -s redis-${REDIS_PV} redis
	echo Please use $(REDIS_SERVER)

$(HIREDIS_LIB): hiredis-${HIREDIS_PV}.tar.gz
	tar -xzvf $<
	make -C hiredis-${HIREDIS_PV}
	ln -s hiredis-${HIREDIS_PV} hiredis

clean:
	rm -rf redis-${REDIS_PV}
	rm -rf hiredis-${HIREDIS_PV}

clean-tar:
		rm -rf *.tar.gz
