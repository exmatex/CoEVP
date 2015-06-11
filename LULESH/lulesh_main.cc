#include "lulesh.h"

#if defined(COEVP_MPI)
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{
#if defined(COEVP_MPI)
   int numRanks ;
   int myRank ;
   MPI_Init(&argc, &argv) ;
   MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
#endif

  Lulesh luleshSystem;

#if defined(COEVP_MPI)
  luleshSystem.go(myRank, numRanks); 
#else
  luleshSystem.go();
#endif

#if defined(COEVP_MPI)
   MPI_Finalize() ;
#endif

  return 0;
}
