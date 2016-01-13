CLEANEXTS   = o a
LOGLIBFILE  = liblogger.a

.PHONY: all
all: $(LOGLIBFILE)

$(LOGLIBFILE): LoggerDB.o Locator.o
	ar ru $@ $^
	ranlib $@

.PHONY: clean
clean:
	for file in $(CLEANEXTS); do rm -f *.$$file; done

LoggerDB.o: LoggerDB.h Logger.h
Locator.o: Locator.h Logger.h
