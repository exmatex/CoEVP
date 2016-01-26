CLEANEXTS   = o a
LOGLIBFILE  = liblogger.a
CXXFLAGS   += -I$(REDIS_LOC) -std=c++11

.PHONY: all
all: $(LOGLIBFILE)

ifeq ($(REDIS), yes)
$(LOGLIBFILE): LoggerDB.o Locator.o
else
$(LOGLIBFILE): Locator.o
endif
	ar ru $@ $^
	ranlib $@

.PHONY: clean
clean:
	for file in $(CLEANEXTS); do rm -f *.$$file; done

ifeq ($(REDIS), yes)
LoggerDB.o: LoggerDB.h Logger.h LoggerDB.cc
endif
Locator.o: Locator.h Logger.h Locator.cc
