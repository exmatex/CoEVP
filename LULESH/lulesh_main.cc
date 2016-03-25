#include <stdexcept>

#if defined(COEVP_MPI)
#include <mpi.h>
#endif

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
  
   int numTaskHandlers = 1;
   int numTasks=2;

   if(std::getenv("NTASKS")!=NULL)
   {
	  numTasks = atoi(std::getenv("NTASKS"));
   }
   if(std::getenv("NHANDLERS")!=NULL)
   {
	  numTaskHandlers = atoi(std::getenv("NHANDLERS"));
   }
#if defined(COEVP_MPI)
   MPI_Init(&argc, &argv) ;
   MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

#if defined(MPI_TASK_POOL)

// I'm really sorry but to avoid initialization headaches, lulesh is both the producer and consumer of work
// so we have to check if we were instantiated by another mpi process
  MPI_Comm mpi_intercomm_parent;
  MPI_Comm_get_parent(&mpi_intercomm_parent);
  MPI_Comm mpi_comm_taskhandler;
  MPI_Comm mpi_intercomm_taskpool;

// create a common intercommunicator between the lulesh domains and the task handlers
  int myDomainID;
  int myHandler;


  if (mpi_intercomm_parent == MPI_COMM_NULL)  
  {

	  int rank, size;
	  MPI_Comm mpi_intercomm_taskhandler;

   	  char **command_argv;
	  
	  command_argv = (char **)malloc(argc+3 * sizeof(char *));
   	  for(int i=0;i<argc;i++)
	  {
		command_argv[i] = argv[i];
	  }
	  command_argv[argc] = (char *)"-E 4";
	  command_argv[argc+1] = (char *)"-H 1";
	  command_argv[argc+2] = NULL;

	  printf("Spawning %d MPI Task Handlers\n", numTaskHandlers);

	  MPI_Comm_spawn((char *)"taskhandler", command_argv, numTaskHandlers, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &mpi_intercomm_taskhandler, MPI_ERRCODES_IGNORE);
 
	  // here it gets complicated. we need to new intracoomunicator including our spawned task handlers, so we can doa collect launch of the lulesh process
	  MPI_Intercomm_merge(mpi_intercomm_taskhandler, 1, &mpi_comm_taskhandler); 

	  MPI_Comm_rank (mpi_comm_taskhandler, &rank);
	  MPI_Comm_size (mpi_comm_taskhandler, &size);
	  printf( "View from Lulesh on intracommunicator  %d of %d\n", rank, size );

	  myDomainID = rank; //this is used to ID me when I request workers

	  // let's tell the task handlers how many tasks we want to spawn

	  MPI_Bcast(&numTasks, 1, MPI_INT, size-1, mpi_comm_taskhandler);

	  // we build a shared intracommunicator, so let's use it to do a collective mpi_spawn on our tasks	

	  printf("LULESH collective spawning %d tasks\n", numTasks);
	  MPI_Comm_spawn(command_argv[0], (command_argv+1), numTasks, MPI_INFO_NULL, size-1, mpi_comm_taskhandler, &mpi_intercomm_taskpool, MPI_ERRCODES_IGNORE);

	  // we have to take part in the collective bcast cool to let all tasks lnow the number of tasks
	  MPI_Bcast(&numTaskHandlers, 1, MPI_INT, MPI_PROC_NULL, mpi_intercomm_taskpool);

	  // I had better figure out my priority task scheduler while I'm at it

	  myHandler = (int) (((float)myRank / (float)numRanks) * (float)numTaskHandlers);

	  printf("Lulesh Rank %d sees that there are %d task handlers. It is affinitised to Task Handler %d\n", myRank, numTaskHandlers, myHandler);

	}
	else
	{


	  // let's broadcast the number of task handlers why not, this is using an intercommunicator so behaves a little difference
  
	  int numTaskHandlers;
 
	  MPI_Bcast(&numTaskHandlers, 1, MPI_INT, 0, mpi_intercomm_parent);
  
	  myHandler = (int) (((float)myRank / (float)numRanks) * (float)numTaskHandlers);
	  printf("Lulesh Task Worker %d sees that there are %d task handlers. It is affinitised to Task Handler %d\n", myRank, numTaskHandlers, myHandler);

		// we need to convince lulesh to ignore mpi domain decomposition, hopefully this hack will do it
 
      numRanks = 1;
	  myRank = 0;

	}

#endif

#endif
  
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
  
  Lulesh luleshSystem;

#if defined(COEVP_MPI)
  #if defined(MPI_TASK_POOL)
  luleshSystem.mpi_comm_taskhandler=mpi_comm_taskhandler;
  luleshSystem.mpi_intercomm_taskpool = mpi_intercomm_taskpool;
  luleshSystem.mpi_intercomm_parent = mpi_intercomm_parent;
  luleshSystem.myDomainID = myDomainID;
  luleshSystem.myHandler = myHandler;
  #endif
#endif

  //  Parse command line optoins
  int  help   = 0;
  
  addArg("help",     'h', 0, 'i',  &(help),                0, "print this message");
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
  addArg("Domain Stop Time", 'D', 1, 'd',  &(domStopTime), 0, "Number of Simulated Seconds to Run For"); 
  processArgs(argc,argv);
  
  if (help) {
    printArgs();
    freeArgs();
    exit(1);
  } 


  // Initialize Taylor cylinder mesh
  luleshSystem.Initialize(myRank, numRanks, edgeElems, heightElems, domStopTime);



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
     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
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
   MPI_Finalize() ;
#endif

  return 0;
}
