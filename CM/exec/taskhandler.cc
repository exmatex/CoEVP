#include <stdio.h>
#include <mpi.h>


int main(int argc, char** argv)
{
  int rank, size;

  MPI_Init (&argc, &argv);	/* starts MPI */

  MPI_Comm mpi_comm_taskhandler;


  // get the parent intercommunicator
  MPI_Comm mpi_intercomm_parent;
  MPI_Comm_get_parent(&mpi_intercomm_parent);
  if (mpi_intercomm_parent != MPI_COMM_NULL)  
  {
    // generate a new intracommunicator from parent

    MPI_Intercomm_merge(mpi_intercomm_parent, 0, &mpi_comm_taskhandler); 

    MPI_Comm_rank (mpi_comm_taskhandler, &rank);
    MPI_Comm_size (mpi_comm_taskhandler, &size);
   
    printf( "Task Handler View from parent communicator  %d of %d\n", rank, size );




  }
  else
  {

    MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm_taskhandler);
  }



  MPI_Comm_rank (mpi_comm_taskhandler, &rank);
  MPI_Comm_size (mpi_comm_taskhandler, &size);
  printf( "Hello world from taskhandler  %d of %d\n", rank, size );


  // first we spawn the taskpool
  MPI_Comm mpi_intercomm_taskpool;
//  MPI_Comm_spawn("/home/vernon/CoEVP/CM/exec/kintask", MPI_ARGV_NULL, 4, MPI_INFO_NULL, 0, mpi_comm_taskhandler, &mpi_comm_taskpool, MPI_ERRCODES_IGNORE);


  MPI_Finalize();
  return 0;
}
