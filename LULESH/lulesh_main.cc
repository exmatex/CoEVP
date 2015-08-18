#include "lulesh.h"

#if defined(COEVP_MPI)
#include <mpi.h>
#endif

//#define USE_ADAPTIVE_SAMPLING

int main(int argc, char *argv[])
{
   int numRanks = 1;
   int myRank = 0;
#if defined(COEVP_MPI)
   MPI_Init(&argc, &argv) ;
   MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
#endif

  Lulesh luleshSystem;

  // Initialize Taylor cylinder mesh
  luleshSystem.Initialize(myRank, numRanks);
  
#ifdef USE_ADAPTIVE_SAMPLING
   bool use_adaptive_sampling = true;
#else
   bool use_adaptive_sampling = false;
#endif

  // Construct fine scale models
  luleshSystem.ConstructFineScaleModel(use_adaptive_sampling);
  
  // Exchange nodal mass
  luleshSystem.ExchangeNodalMass();

  // Simulate 
  luleshSystem.go(myRank, numRanks, use_adaptive_sampling); 

#if defined(COEVP_MPI)
   MPI_Finalize() ;
#endif

  return 0;
}
