FLANN_PV=1.8.1
FLANN_DIR=flann-${FLANN_PV}-src
FLANN_LIB=$(FLANN_DIR)/lib/libflann_cpp.so

all: ${FLANN_LIB}

$(FLANN_DIR).tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://github.com/mariusmuja/flann/archive/${FLANN_PV}-src.tar.gz
else
	curl -o $@ -L https://github.com/mariusmuja/flann/archive/${FLANN_PV}-src.tar.gz
endif

$(FLANN_LIB): ${FLANN_DIR}.tar.gz
	tar -xzvf $<
	cd ${FLANN_DIR} && cmake -DPDFLATEX_COMPILER=NOTFOUND .
	make -C ${FLANN_DIR}
	ln -s ${FLANN_DIR} flann
	echo Please use $(FLANN_LIB)

clean:
	rm -rf ${FLANN_DIR} flann

clean-tar:
		rm -rf *.tar.gz
