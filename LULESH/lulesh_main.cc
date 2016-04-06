#include <stdexcept>

#if defined(COEVP_MPI)
#include <mpi.h>
#endif

#include <chrono>
#include "lulesh.h"
#include "SingletonDB.h"
#include "ModelDB_SingletonDB.h"

#if defined(LOGGER)      // CoEVP Makefile enforces assert: LOGGER=REDIS=yes
#include "LoggerDB.h"    // Includes Logger base class too
#include "Locator.h"
#endif


//  Command line option parsing (using Sriram code from old days)
#include "cmdLineParser.h"

int main(int argc, char *argv[])
{
   int numRanks = 1;
   int myRank = 0;

  
   int numTaskHandlers = 0;
   int numTasks=0;
   int numCoEVP=0;

   int myDomainID = 0;
   
#if defined(COEVP_MPI)
   MPI_Init(&argc, &argv) ;
   MPI_Comm myComm;
   MPI_Comm_dup(MPI_COMM_WORLD, &myComm);
   
   MPI_Comm_size(myComm, &numRanks) ;
   MPI_Comm_rank(myComm, &myRank) ;


   if(std::getenv("NTASKS")!=NULL)
   {
      numTasks = atoi(std::getenv("NTASKS"));
      if(std::getenv("NHANDLERS")!=NULL)
      {
      	  numTaskHandlers = atoi(std::getenv("NHANDLERS"));
		  numCoEVP = numRanks - numTaskHandlers - numTasks;
		  if(myRank==0)
		  {
		  	std::cout << "Task Pool: CoEVP Ranks: " << numCoEVP << " TaskHandlers: " << numTaskHandlers << " Task Workers: " << numTasks << std::endl;
		  }
		  
		  if(numCoEVP<1 || numTaskHandlers<1 || numTasks<1)
		  {
				std::cout << "This clearly isn't going to work is it." << std::endl;
				exit(1);
		  }

      }
	  
   }




// to avoid initialization headaches, lulesh is both the producer and consumer of work
// so we have to check if we were instantiated by another mpi process
  MPI_Comm mpi_intercomm_taskhandler;
  MPI_Comm mpi_comm_taskhandler;
  MPI_Comm mpi_intercomm_taskpool;
  MPI_Comm mpi_comm_taskpool;
  MPI_Comm CoEVP_comm;

// create a common intercommunicator between the lulesh domains and the task handlers
  int myHandler = 0;


  if (numTaskHandlers)  
  {


	// this is the refactor for the non comm_spawn implementation
	// order in ranks is handlers, coevp, task workers
	// workers and coevp will see this code, not the handlers
	// let's build our communicators accordingly

	  MPI_Group world_group;
	  MPI_Comm_group(myComm, &world_group);

	
      MPI_Group task_group;
	  MPI_Group taskhandler_group;
      MPI_Group CoEVP_group;
	
	  int task_range[1][3];
	  int taskhandler_range[1][3];
	  int CoEVP_range[1][3];
      
	  task_range[0][1] = numRanks - 1;
	  task_range[0][0] = numRanks - numTasks;
	  task_range[0][2] = 1;

	  taskhandler_range[0][1] = task_range[0][0] - 1;
	  taskhandler_range[0][0] = 0;
	  taskhandler_range[0][2] = 1;
     
      CoEVP_range[0][0] = numTaskHandlers;
      CoEVP_range[0][1] = numTaskHandlers+numCoEVP-1;
      CoEVP_range[0][2] = 1;
     
	  MPI_Group_range_incl(world_group, 1, task_range, &task_group);
	  MPI_Group_range_incl(world_group, 1, taskhandler_range, &taskhandler_group);
	  MPI_Group_range_incl(world_group, 1, CoEVP_range, &CoEVP_group);

	  MPI_Comm_create_group(myComm, task_group, 0, &mpi_comm_taskpool);
	  MPI_Comm_create_group(myComm, taskhandler_group, 0, &mpi_comm_taskhandler);
	  MPI_Comm_create_group(myComm, CoEVP_group, 0, &CoEVP_comm);
	
      // we built two communicators one for task workers and one for everything else
	  // now let's build the intercommunicators between the groups

	  std::cout << "Tasks: " << task_range[0][0] << ".." << task_range[0][1] << std::endl;
	  std::cout << "Handler and CoEVP: " << taskhandler_range[0][0] << ".." << taskhandler_range[0][1] << std::endl;
	  std::cout << "CoEVP: " << CoEVP_range[0][0] << ".." << CoEVP_range[0][1] << std::endl;      
	  

	  if(myRank > task_range[0][0]-1)
	  {
		// i am a task worker
		myHandler = (int) (((float) (myRank - task_range[0][0]) / (float)numTasks) * (float)numTaskHandlers);

		std::cout << "I am task worker: " << myRank - task_range[0][0] << ", with global rank: " << myRank << std::endl;
		MPI_Comm_size(mpi_comm_taskpool, &numRanks);
		MPI_Comm_rank(mpi_comm_taskpool, &myRank);
		std::cout << "I am task worker: " << myRank << ", of size: " << numRanks << std::endl;


		MPI_Intercomm_create(mpi_comm_taskpool, numRanks-1, myComm, 0, 1, &mpi_intercomm_taskhandler);
                 
		// don't want the tasks doing any confused collectives
//		MPI_Comm_dup(MPI_COMM_SELF, &myComm);
		myComm = MPI_COMM_SELF;

		printf("Lulesh Task Worker %d sees that there are %d task handlers. It is affinitised to Task Handler %d\n", myRank, numTaskHandlers, myHandler);

		MPI_Comm_size(myComm, &numRanks);
		MPI_Comm_rank(myComm, &myRank);


	  }
	  else
	  {

		// i am a coevp rank or taskhandler
        std::cout << "I am taskhandler (or coevp rank): " << myRank << std::endl;
		MPI_Intercomm_create(mpi_comm_taskhandler, 0, myComm, numRanks-1, 1, &mpi_intercomm_taskpool);
      
		
        mpi_comm_taskpool = MPI_COMM_NULL;
        mpi_intercomm_taskhandler = MPI_COMM_NULL;
		// if I am a regular coevp rank we need to clean up the communicator so coevp doesn't get confused
		if(myRank>numTaskHandlers-1)
		{


	
            myComm = CoEVP_comm;
			myDomainID = myRank; // this is my global rank 

			MPI_Comm_size(myComm, &numRanks);
			MPI_Comm_rank(myComm, &myRank);

	  		//myHandler = (int) (((float)myRank / (float)numRanks) * (float)numTaskHandlers);
		    myHandler = (int) (((float) (myRank) / (float)(numRanks)) * (float)numTaskHandlers);            

			printf("CoEVP Rank %d sees that there are %d ranks and %d task handlers. It is affinitised to Task Handler %d\n", myRank, numRanks, numTaskHandlers, myHandler);


		}
		else
		{
            
			// task handler logic goes here
            std::cout << "Starting CoEVP Task Handler" << std::endl;
	        int lulesh_work_probe;
    	    int task_worker_probe;
        	MPI_Status mpi_status;
		    MPI_Request mpi_request;
    	    std::list<int> lulesh_work;
        	std::list<int> task_worker;

	        while(1)
    	    {

                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm_taskhandler, &lulesh_work_probe, &mpi_status);

                if(lulesh_work_probe)
                {
                        // there is some incoming work
    			        int lulesh_work_id;
                        MPI_Recv(&lulesh_work_id, 1, MPI_INT, MPI_ANY_SOURCE, 1, mpi_comm_taskhandler, &mpi_status);
                        //std::cout << "TaskHandler received work: " << lulesh_work_id << std::endl;
                        // do we have any registered idle workers?
                        if(task_worker.size()>0)
                        {
                                // pair up the task_worker with the work
                                int task_worker_id = task_worker.front();
                                task_worker.pop_front();
                		        MPI_Send(&lulesh_work_id, 1, MPI_INT, task_worker_id, 3, mpi_intercomm_taskpool);//, &mpi_request);
                        }
                        else if(numTaskHandlers>1)
                        {
                                // we received work but we don't have any workers to assign it to so let's round robin to other taskhandlers (if they exist)
//                              std::cout<<"Load balancing task"<<std::endl;
                                if(myRank == numTaskHandlers-1)
                                {
                                        //we're the last rank, wrap around to handler 0
                                        MPI_Send(&lulesh_work_id, 1, MPI_INT, 0, 1, mpi_comm_taskhandler);//, &mpi_request);
                                }
                                else
                                {
                                        //send task to neighbouring handler
                                        MPI_Send(&lulesh_work_id, 1, MPI_INT, myRank+1, 1, mpi_comm_taskhandler);//, &mpi_request);
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
            			int task_worker_id;
                        MPI_Recv(&task_worker_id, 1, MPI_INT, MPI_ANY_SOURCE, 2, mpi_intercomm_taskpool, &mpi_status);
                        if(lulesh_work.size()>0)
                        {
                                //imediately pair up with idle task_worker
                                int lulesh_work_id = lulesh_work.front();
                                lulesh_work.pop_front();
                                MPI_Send(&lulesh_work_id, 1, MPI_INT, task_worker_id, 3, mpi_intercomm_taskpool);//, &mpi_request);
                        }
                        else
                        {
                                task_worker.push_back(task_worker_id);
//                              std::cout << "Number of queued task workers:" << task_worker.size() << std::endl;       
                        }
                }

	        }
		}

	  }
  }
  else
  {
	// just a regular mpi run
	mpi_intercomm_taskhandler = MPI_COMM_NULL;
  }

#endif
  int  timer = 0;
  int  sampling = 0;              //  By default, use adaptive sampling (but compiled in)
  int  redising = 0;              //  By default, do not use REDIS for database
  int  posixing = 0;              // POSIX Key/Value store
  int  hioing = 0;                //  HIO global store
  int  global_ns = 0;             //  By default, do not use a global earest neighbor
  int  flanning = 0;              //  By default, do not use FLANN for nearest neighbor search
  int  logging = 0;               //  By default, do not use FLANN for nearest neighbor search
  int  flann_n_trees = 1;         // Default can be overridden using command line
  int  flann_n_checks = 20;       // Default can be overridden using command line
  int  file_parts = 0;
  int  debug_topology = 0;
  int  visit_data_interval = 0;  // Set this to 0 to disable VisIt data writing
  int  distributed_redis = 0;
  char logdb[1024] = {0};        // host and port of logging databse (e.g. cn1:6379)
  int heightElems = 26;
  int edgeElems = 16;
  double domStopTime = 1.e-1;
  int simStopCycle = 0;
  
  Lulesh luleshSystem;

  luleshSystem.myDomainID = myDomainID;

#if defined(COEVP_MPI)

  luleshSystem.mpi_comm_taskhandler=mpi_comm_taskhandler;
  luleshSystem.mpi_comm_taskpool=mpi_comm_taskpool;
  luleshSystem.mpi_intercomm_taskpool = mpi_intercomm_taskpool;
  luleshSystem.mpi_intercomm_taskhandler = mpi_intercomm_taskhandler;
  luleshSystem.myHandler = myHandler;
  luleshSystem.myComm = myComm;

#endif

  //  Parse command line optoins
  int  help   = 0;
  
  addArg("help",     'h', 0, 'i',  &(help),                0, "print this message");
  addArg("timer",    'a', 1, 'i',  &(timer),               0, "use timing and write to output file");
  addArg("sample",   's', 0, 'i',  &(sampling),            0, "use adaptive sampling");
  addArg("redis",    'r', 0, 'i',  &(redising),            0, "use REDIS library");
  addArg("posix",    'x', 0, 'i',  &(posixing),            0, "use POSIX library");
  addArg("hio",      'b', 0, 'i',  &(hioing),              0, "use HIO library");
  addArg("globalns" ,'g', 0, 'i',  &(global_ns),           0, "use global neighbor search/data store");
  addArg("flann",    'f', 0, 'i',  &(flanning),            0, "use FLANN library");
  addArg("n_trees",  't', 1, 'i',  &(flann_n_trees),       0, "number of FLANN trees");
  addArg("n_checks", 'c', 1, 'i',  &(flann_n_checks),      0, "number of FLANN checks");
  addArg("parts",    'p', 1, 'i',  &(file_parts),          0, "number of file parts");
  addArg("visitint", 'v', 1, 'i',  &(visit_data_interval), 0, "visit output interval");
  addArg("debug",    'd', 0, 'i',  &(debug_topology),      0, "add debug info to SILO");
  addArg("distributed_redis", 'R', 0, 'i', &(distributed_redis), 0, "use distributed REDIS via twemproxy");
  addArg("log",      'l', 1, 's',  &(logdb),   sizeof(logdb), "log to REDIS at hostname:port");
  addArg("Height Elems", 'H', 1, 'i', &(heightElems), 0, "Number of height elements to solve for");
  addArg("Edge Elems", 'E', 1, 'i', &(edgeElems), 0, "Number of height elements to solve for");
  addArg("Domain Stop Cycle", 'C', 1, 'i', &(simStopCycle), 0, "Number of Simulated Cycles to Run For");
  addArg("Domain Stop Time", 'D', 1, 'd',  &(domStopTime), 0, "Number of Simulated Seconds to Run For"); 
  processArgs(argc,argv);
  
  if (help) {
    printArgs();
    freeArgs();
    exit(1);
  } 


  // Initialize Taylor cylinder mesh
  luleshSystem.Initialize(myRank, numRanks, edgeElems, heightElems, domStopTime, simStopCycle, timer);



  if (sampling) {
    printf("Using adaptive sampling...\n");
  } else {
    if (redising||distributed_redis||flanning||global_ns) {
      throw std::runtime_error("--redis/--distributed_redis/--flann/--globalns needs --sample");
    }
  }
  if (redising) {
    printf("Using Redis library...\n");
    if(distributed_redis)
      printf("Using Distributed Redis (twemproxy)...\n");
  }
  if (posixing) 
  {
    printf("Using POSIX library...\n");
  }
  if (hioing) 
  {
    printf("Using HIO library...\n");
  }
  if (flanning) {
    printf("Using FLANN library...\n");
    printf("   flann_n_trees: %d\n", flann_n_trees);
    printf("   flann_n_checks: %d\n", flann_n_checks);
  }
  if (visit_data_interval != 0){
#ifndef SILO
      throw std::runtime_error("--redis/--flann/--globalns needs --sample"); 
#endif
  }
  if (strlen(logdb) != 0) {
    logging = 1;
  }
  freeArgs();

   /*************************************/
   /* Initialize ModelDB Interface      */
   /*************************************/
   ModelDatabase * global_modelDB = nullptr;
   ApproxNearestNeighbors* global_ann = nullptr;
   if(sampling)
   {
      if(redising){
        if(distributed_redis)
          SingletonDB::getInstance(SingletonDBBackendEnum::DIST_REDIS_DB);
        else
          SingletonDB::getInstance(SingletonDBBackendEnum::REDIS_DB);
        global_modelDB = new ModelDB_SingletonDB();
      }
      // Add in new (global) backends that require initialization here
      // Optionally put a !variable in the else if
	  else if(posixing) {
        SingletonDB::getInstance(SingletonDBBackendEnum::POSIX_DB);
        global_modelDB = new ModelDB_SingletonDB();
	  }
	  else if(hioing) {
        SingletonDB::getInstance(SingletonDBBackendEnum::HIO_DB);
        global_modelDB = new ModelDB_SingletonDB();
	  }
      else if(!redising && global_ns){
        SingletonDB::getInstance(SingletonDBBackendEnum::HASHMAP_DB);
        global_modelDB = new ModelDB_SingletonDB();
      }
      else
      {
        //Do nothing because we don't have to initialize a global_modelDB
      }
   }

#if defined(LOGGER)
   // Initialize logging to REDIS database
   Locator::initialize();             // make sure a dummy logger exists
   if (logging) {
#if defined(COEVP_MPI)
     int my_rank;
     MPI_Comm_rank(myComm, &my_rank);
     char my_node[MPI_MAX_PROCESSOR_NAME];
     int name_len;
     MPI_Get_processor_name(my_node, &name_len);
     LoggerDB  *logger_db = new LoggerDB(logdb, std::string(my_node), my_rank);
#else
     LoggerDB  *logger_db = new LoggerDB(logdb);
#endif
     Locator::provide(logger_db);
   }

   Logger  &logger = Locator::getLogger();
#endif

  // Construct fine scale models
  luleshSystem.ConstructFineScaleModel(sampling,global_modelDB,global_ann,flanning,flann_n_trees,flann_n_checks,global_ns);
  
  // Exchange nodal mass
  luleshSystem.ExchangeNodalMass();

  // Simulate 
  luleshSystem.go(myRank,numRanks,sampling,visit_data_interval,file_parts,debug_topology);

  // Only do this is we have actually opened a REDIS connection.
#if defined(LOGGER)
  if (logging) {
    delete(&logger); 
  }
#endif
  
#if defined(COEVP_MPI)
//   MPI_Finalize() ;
#endif

  return 0;
}
