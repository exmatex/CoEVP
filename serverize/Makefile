PROTOBUF_PV=2.6.1
PROTOBUF_DIR=protobuf-${PROTOBUF_PV}
PROTOBUF_INSTALL=protobuf
PROTOBUF_BIN=${PROTOBUF_INSTALL}/bin/protoc

CIRCLE_PV=0.2.1-rc.1
CIRCLE_DIR=libcircle-${CIRCLE_PV}
CIRCLE_INSTALL=circle
CIRCLE_LIB=${CIRCLE_INSTALL}/lib/libcircle.a

all: ${PROTOBUF_BIN} ${CIRCLE_LIB}

protobuf: ${PROTOBUF_BIN}

$(PROTOBUF_DIR).tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://github.com/google/protobuf/releases/download/v${PROTOBUF_PV}/protobuf-${PROTOBUF_PV}.tar.gz
else
	curl -o $@ -L https://github.com/google/protobuf/releases/download/v${PROTOBUF_PV}/protobuf-${PROTOBUF_PV}.tar.gz
endif

$(PROTOBUF_BIN): ${PROTOBUF_DIR}.tar.gz
	tar -xzvf $<
	mkdir -p ${PROTOBUF_INSTALL}
	cd ${PROTOBUF_DIR}; ./configure --prefix=${PWD}/${PROTOBUF_INSTALL}; make; make install

$(CIRCLE_DIR).tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://github.com/hpc/libcircle/releases/download/${CIRCLE_PV}/libcircle-${CIRCLE_PV}.tar.gz
else
	curl -o $@ -L https://github.com/hpc/libcircle/releases/download/${CIRCLE_PV}/libcircle-${CIRCLE_PV}.tar.gz
endif

$(CIRCLE_LIB): ${CIRCLE_DIR}.tar.gz
	tar -xzvf $<
	cd ${CIRCLE_DIR}; ./configure --prefix=${PWD}/${CIRCLE_INSTALL}; make; make install

clean:
	rm -rf ${PROTOBUF_DIR} ${CIRCLE_DIR} protobuf circle

clean-tar:
		rm -rf *.tar.gz
