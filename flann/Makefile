FLANN_PV=1.8.1

FLANN_LIB=build/lib/libflann_cpp.so

all: ${FLANN_LIB}

flann-${FLANN_PV}-src.tar.gz:
	wget -O $@ https://github.com/mariusmuja/flann/archive/${FLANN_PV}-src.tar.gz

$(FLANN_LIB): flann-${FLANN_PV}-src.tar.gz
	tar -xzvf $<
	mkdir build
	cd build && cmake ../flann-${FLANN_PV}-src
	make -C build
	ln -s flann-${FLANN_PV}-src flann
	echo Please use $(FLANN_LIB)

clean:
	rm -rf flann-${FLANN_PV}-src flann build

clean-tar:
		rm -rf *.tar.gz
