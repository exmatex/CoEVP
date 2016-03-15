#include <stdio.h>
#include <mpi.h>


int main(int argc, char** argv)
{
  int rank, size, numTasks;

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

    MPI_Bcast(&numTasks, 1, MPI_INT, size-1, mpi_comm_taskhandler);

    printf("Task Handler collective launching %d tasks\n", numTasks);

    // now we do a collective call to initialize the task pool (along with lulesh so the tasks can do call backs to lulesh when done)
    MPI_Comm mpi_intercomm_taskpool;
    MPI_Comm_spawn("/home/vernon/CoEVP/CM/exec/kintask", MPI_ARGV_NULL, numTasks, MPI_INFO_NULL, size-1, mpi_comm_taskhandler, &mpi_intercomm_taskpool, MPI_ERRCODES_IGNORE);


  }
  else
  {

    MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm_taskhandler);
  }



  MPI_Comm_rank (mpi_comm_taskhandler, &rank);
  MPI_Comm_size (mpi_comm_taskhandler, &size);





  MPI_Finalize();
  return 0;
}
