#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <libcircle.h>
#include <mpi.h>

int timestep = 0;


//  Get functions parameters (unpack protobuf-encoded string), call
//  function, and write results to database.
void processTasks(CIRCLE_handle *handle) {
  char taskString[96];
  handle->dequeue(&taskString[0]);

  int numRanks;
  int myRank;
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
  char my_node[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(my_node, &name_len);
  std::cout << my_node << "/" << myRank << "/" << numRanks << " task string: " << taskString << std::endl;
}


//  Setup (via protobuf) parameters to function call.
void buildTasks(CIRCLE_handle *handle) {
  char taskString[96];
  for(int i=0; i<5; i++) {
    sprintf(taskString, "timestep: %d,  task: %d", timestep, i);
    handle->enqueue(taskString);
  }
}


void doCircleTasks() {
  CIRCLE_cb_create(&buildTasks);
  CIRCLE_begin();
}


int main(int argc, char ** argv)
{
  int rank = CIRCLE_init(argc, argv, CIRCLE_DEFAULT_FLAGS);
  CIRCLE_cb_process(&processTasks);
  CIRCLE_enable_logging(CIRCLE_LOG_ERR);

  if(rank != 0) {
    for(int i=0; i<3; i++) {
      CIRCLE_begin();
      //  MPI_Barrier(MPI_COMM_WORLD);
      //  At this point, all libcircle tasks are finished(?)
    }
    CIRCLE_finalize();
    return 0;
  }

  for(int i=0; i<3; i++) {
    timestep = i;
    doCircleTasks();
    //  Is this even necessary, or does CIRCLE_begin() block here too?
    //  MPI_Barrier(MPI_COMM_WORLD);
    //  Read results from database until it's "empty". Remember that it is
    //  eventually consistent so there is a chance that we will not get all
    //  results. There is a way to check it (using atomic INCR to count
    //  executions), but no way to ensure that we have all results. Ugh.
    sleep(5);
    std:: cout << "------" << std::endl;
  }
  CIRCLE_finalize();
  return 0;
}
