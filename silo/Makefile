SILO_PV=4.10.2
SILO_LIB=silo/lib/libsiloh5.a

all: $(SILO_LIB)

silo-${SILO_PV}.tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-$(SILO_PV)/silo-$(SILO_PV).tar.gz 
else
	curl -o $@ -L https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-$(SILO_PV)/silo-$(SILO_PV).tar.gz 
endif

$(SILO_LIB): silo-${SILO_PV}.tar.gz
	tar -xzvf $<
	cd silo-${SILO_PV} && ./configure --prefix=$(CURDIR)/silo --with-hdf5=/usr/include,/usr/lib --enable-shared
	make -C silo-${SILO_PV}
	make -C silo-${SILO_PV} install

clean:
	rm -rf silo-${SILO_PV} silo

clean-tar:
		rm -rf *.tar.gz
