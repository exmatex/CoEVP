ifneq ($(REDIS), yes)
$(info ** ERROR: logging requires REDIS, set REDIS=yes in top Makefile **)
endif

CLEANEXTS   = o a
LOGLIBFILE  = liblogger.a
CXXFLAGS   += -I$(REDIS_LOC) -std=c++11

.PHONY: all
all: $(LOGLIBFILE)

$(LOGLIBFILE): LoggerDB.o Locator.o
	ar ru $@ $^
	ranlib $@

.PHONY: clean
clean:
	for file in $(CLEANEXTS); do rm -f *.$$file; done

LoggerDB.o: LoggerDB.h Logger.h LoggerDB.cc
Locator.o: Locator.h Logger.h Locator.cc
