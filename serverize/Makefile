PROTOBUF_PV=2.6.1
PROTOBUF_DIR=protobuf-${PROTOBUF_PV}
PROTOBUF_INSTALL=protobuf
PROTOBUF_BIN=${PROTOBUF_INSTALL}/bin/protoc

all: ${PROTOBUF_BIN}

$(PROTOBUF_DIR).tar.gz:
ifeq (, $(shell which curl))
	wget -O $@ https://github.com/google/protobuf/releases/download/v${PROTOBUF_PV}/protobuf-${PROTOBUF_PV}.tar.gz
else
	curl -o $@ -L https://github.com/google/protobuf/releases/download/v${PROTOBUF_PV}/protobuf-${PROTOBUF_PV}.tar.gz
endif

$(PROTOBUF_BIN): ${PROTOBUF_DIR}.tar.gz
	tar -xzvf $<
	mkdir ${PROTOBUF_INSTALL}
	cd ${PROTOBUF_DIR}; ./configure --prefix=${PWD}/${PROTOBUF_INSTALL}; make; make install

clean:
	rm -rf ${PROTOBUF_DIR} protobuf

clean-tar:
		rm -rf *.tar.gz
