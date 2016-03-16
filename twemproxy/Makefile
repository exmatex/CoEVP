REDIS_PV=2.8.22
HIREDIS_PV=0.13.3

TWEMPROXY_SERVER=twemproxy-master/src/nutcracker

REDIS_SERVER=redis-${REDIS_PV}/src/redis-server
HIREDIS_LIB=hiredis-${HIREDIS_PV}/libhiredis.a

all: $(TWEMPROXY_SERVER)

master.tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://github.com/twitter/twemproxy/archive/master.tar.gz $(WGETFLAG)
else
	curl -o $@ -L https://github.com/twitter/twemproxy/archive/master.tar.gz $(CURLFLAG)
endif

$(TWEMPROXY_SERVER): master.tar.gz
	tar -xzvf $<
	cd twemproxy-master; autoreconf -fvi; ./configure
	${MAKE} -C twemproxy-master

clean:
	${MAKE} -C twemproxy-master clean
	${MAKE} -C twemproxy-master distclean

clean-tar:
		rm -rf *.tar.gz
