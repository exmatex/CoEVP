#ifndef __SILODUMP_H__
#define __SILODUMP_H__

#include "lulesh.h"

void DumpDomain(LuleshNS::Domain *domain, int myRank, int numProcs, int fileParts, int sampling, int debug_topology);

#endif
