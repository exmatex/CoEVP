FLANN_PV=1.8.1
FLANN_DIR=flann-${FLANN_PV}-src
FLANN_HEADER=$(FLANN_DIR)/src/cpp/flann/flann.hpp

all: ${FLANN_HEADER}

$(FLANN_DIR).tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://github.com/mariusmuja/flann/archive/${FLANN_PV}-src.tar.gz
else
	curl -o $@ -L https://github.com/mariusmuja/flann/archive/${FLANN_PV}-src.tar.gz
endif

$(FLANN_HEADER): ${FLANN_DIR}.tar.gz
	tar -xzvf $<
	ln -sf ${FLANN_DIR} flann
	touch $(FLANN_HEADER)
	echo Please use $(FLANN_HEADER)

clean:
	rm -rf ${FLANN_DIR} flann

clean-tar:
		rm -rf *.tar.gz
