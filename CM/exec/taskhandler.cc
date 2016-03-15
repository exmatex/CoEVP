#include <stdio.h>
#include <mpi.h>


int main(int argc, char** argv)
{
  int rank, size;

  MPI_Init (&argc, &argv);	/* starts MPI */

  MPI_Comm mpi_comm_taskhandler;
  MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm_taskhandler);

  MPI_Comm_rank (mpi_comm_taskhandler, &rank);
  MPI_Comm_size (mpi_comm_taskhandler, &size);
  printf( "Hello world from process %d of %d\n", rank, size );

  MPI_Comm mpi_comm_taskpool;
  MPI_Comm_spawn("kintask", MPI_ARGV_NULL, 4, MPI_INFO_NULL, 0, mpi_comm_taskhandler, &mpi_comm_taskpool, MPI_ERRCODES_IGNORE);



  MPI_Finalize();
  return 0;
}
