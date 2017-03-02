SILO_PV=4.10.2
SILO_LIB=silo/lib/libsiloh5.a

# Conditionally set defaults for --with-hdf5 if COEVP_SILO_WITH_HDF5 is not already set.
COEVP_SILO_WITH_HDF5 ?= "/usr/include,/usr/lib"

all: $(SILO_LIB)

silo-${SILO_PV}.tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-$(SILO_PV)/silo-$(SILO_PV).tar.gz $(WGETFLAG)
else
	curl -o $@ -L https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-$(SILO_PV)/silo-$(SILO_PV).tar.gz $(CURLFLAG)
endif

$(SILO_LIB): silo-${SILO_PV}.tar.gz
	tar -xzvf $<
	cd silo-${SILO_PV} && ./configure --prefix=$(CURDIR)/silo --with-hdf5=$(COEVP_SILO_WITH_HDF5) --enable-shared --disable-fortran
	make -C silo-${SILO_PV}
	make -C silo-${SILO_PV} install

clean:
	rm -rf silo-${SILO_PV} silo

clean-tar:
		rm -rf *.tar.gz
