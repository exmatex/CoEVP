# makefile for tensor classes

# switches for gnu g++
CXX=	g++
#for g++ 2.9x define NO_STD_LIMITS
#CXXFLAGS=	-O2 -DNO_STD_LIMITS -Wall
CXXFLAGS=	-O2 -Wall
LDXFLAGS=	-O2
CFLAGS=	-O2
LDFLAGS=	-O2
CC=	gcc

# switches for Sun CC
#CXX=	CC
#CXXFLAGS=	-fast
#LDXFLAGS=	-fast
#CFLAGS=	-fast
#LDFLAGS=	-fast
#CC=	cc

# switches for HP aCC
#CXX=   aCC
#CXXFLAGS=  -fast -AA
#LDXFLAGS=  -fast -AA
#CFLAGS=	-fast
#LDFLAGS=   -fast
#CC=    cc

.SUFFIXES: .cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c $*.cc

OBJS=	tensor.o xtensor.o tutils.o basetest.o xtest.o

default:	tensor.o

all:	basetest xtest

basetest:	tensor.o basetest.o
	$(CXX) $(LDXFLAGS) -o basetest tensor.o basetest.o -lm
	
xtest:	tensor.o xtensor.o tutils.o xtest.o
	$(CXX) $(LDXFLAGS) -o xtest tensor.o xtensor.o tutils.o xtest.o -lm
	
test: basetest xtest testdiff
	-./basetest > basetest.out
	-./testdiff basetest.ref basetest.out
	-./xtest > xtest.out
	-./testdiff xtest.ref xtest.out

testdiff:	testdiff.o
	$(CC) $(LDFLAGS) -o testdiff testdiff.o -lm

clean:	
	/bin/rm -f $(OBJS) basetest xtest basetest.out xtest.out testdiff.o testdiff

depend:
	makedepend -- $(CXXFLAGS) -- $(OBJS:o=cc)

# DO NOT DELETE

tensor.o: tensor.h
xtensor.o: xtensor.h tensor.h tutils.h
tutils.o: tutils.h
basetest.o: tensor.h
xtest.o: xtensor.h tensor.h
