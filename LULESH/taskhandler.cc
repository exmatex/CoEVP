#include <stdio.h>
#include <mpi.h>
#include <list>  

int main(int argc, char** argv)
{
  int localrank, numTaskHandlers, rank, size, numTasks;

  MPI_Init (&argc, &argv);	/* starts MPI */
 
  MPI_Comm_rank (MPI_COMM_WORLD, &localrank);
  MPI_Comm_size (MPI_COMM_WORLD, &numTaskHandlers);

  MPI_Comm mpi_comm_taskhandler;


  printf("Total task handlers: %d\n", numTaskHandlers);

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

       char **command_argv;
        
        command_argv = (char **)malloc(argc+3 * sizeof(char *));
        for(int i=0;i<argc;i++)
        {
            std::cout << argv[i] << std::endl;
            command_argv[i] = argv[i];
        }
		command_argv[argc] = "-E 4";
        command_argv[argc+1] = "-H 1";
        command_argv[argc+2] = NULL;


//    MPI_Comm_spawn("/home/vernon/CoEVP/LULESH/lulesh", command_argv, numTasks, MPI_INFO_NULL, size-1, mpi_comm_taskhandler, &mpi_intercomm_taskpool, MPI_ERRCODES_IGNORE);
    MPI_Comm_spawn(command_argv[1], command_argv+2, numTasks, MPI_INFO_NULL, size-1, mpi_comm_taskhandler, &mpi_intercomm_taskpool, MPI_ERRCODES_IGNORE);
//      MPI_Comm_spawn("/home/vernon/CoEVP/CM/exec/kintask", MPI_ARGV_NULL, numTasks, MPI_INFO_NULL, size-1, mpi_comm_taskhandler, &mpi_intercomm_taskpool, MPI_ERRCODES_IGNORE);

    // collective broadcast of number of task handlers to all tasks, MPI_ROOT as we are using an intercommunicator

    if(rank==0)
    {
        MPI_Bcast(&numTaskHandlers, 1, MPI_INT, MPI_ROOT, mpi_intercomm_taskpool);
    }
    else
    {
	    MPI_Bcast(&numTaskHandlers, 1, MPI_INT, MPI_PROC_NULL, mpi_intercomm_taskpool);
    }

    // at this point we are ready to start pairing tasks with lulesh ranks


    int lulesh_work_id;
    int task_worker_id;

	int lulesh_work_probe;
	int task_worker_probe;
	MPI_Status mpi_status;
	std::list<int> lulesh_work;
	std::list<int> task_worker;

	while(1)
	{

		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm_taskhandler, &lulesh_work_probe, &mpi_status);	
		if(lulesh_work_probe)
		{
			// there is some incoming work
			MPI_Recv(&lulesh_work_id, 1, MPI_INT, MPI_ANY_SOURCE, 1, mpi_comm_taskhandler, &mpi_status);
			// do we have any registered idle workers?
			if(task_worker.size()>0)
			{
				// pair up the task_worker with the work
				task_worker_id = task_worker.front();
				task_worker.pop_front();
		        MPI_Send(&lulesh_work_id, 1, MPI_INT, task_worker_id, 3, mpi_intercomm_taskpool);
			}
			else if(numTaskHandlers>1)
			{
				// we received work but we don't have any workers to assign it to so let's round robin to other taskhandlers (if they exist)
				if(localrank == numTaskHandlers-1)
				{
					//we're the last rank, wrap around to handler 0
					MPI_Send(&lulesh_work_id, 1, MPI_INT, 0, 1, mpi_comm_taskhandler);
				}
				else
				{
					//send task to neighbouring handler
					MPI_Send(&lulesh_work_id, 1, MPI_INT, localrank+1, 1, mpi_comm_taskhandler);
				}
			}
			else
			{
				// we're on our own, so let's queue the work
				lulesh_work.push_back(lulesh_work_id);
			}
		}
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG,  mpi_intercomm_taskpool, &task_worker_probe, &mpi_status);	
		if(task_worker_probe)
		{
			MPI_Recv(&task_worker_id, 1, MPI_INT, MPI_ANY_SOURCE, 2, mpi_intercomm_taskpool, &mpi_status);
			if(lulesh_work.size()>0)
			{
				//imediately pair up with idle task_worker
				lulesh_work_id = lulesh_work.front();
				lulesh_work.pop_front();
				MPI_Send(&lulesh_work_id, 1, MPI_INT, task_worker_id, 3, mpi_intercomm_taskpool);
			}
			else
			{
				task_worker.push_back(task_worker_id);
//				std::cout << "Number of queued task workers:" << task_worker.size() << std::endl;	
			}
		}

	}

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
