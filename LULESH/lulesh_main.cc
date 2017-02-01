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


#include "legion.h"

using namespace Legion;

enum TaskID
{
  TOP_LEVEL_TASK_ID,
  MPI_INTEROP_TASK_ID,
  WORKER_TASK_ID,
};

// Here is our global MPI-Legion handshake
// You can have as many of these as you 
// want but the common case is just to
// have one per Legion-MPI rank pair
MPILegionHandshake handshake;

// Have a global static number of iterations for
// this example, but you can easily configure it
// from command line arguments which get passed 
// to both MPI and Legion
const int total_iterations = 10;

void worker_task(const Task *task, 
                 const std::vector<PhysicalRegion> &regions,
                 Context ctx, Runtime *runtime)
{
  printf("Legion Doing Work in Rank %lld\n", 
          task->parent_task->index_point[0]);
}

void mpi_interop_task(const Task *task, 
                      const std::vector<PhysicalRegion> &regions,
                      Context ctx, Runtime *runtime)
{
  printf("Hello from Legion MPI-Interop Task %lld\n", task->index_point[0]);
  for (int i = 0; i < total_iterations; i++)
  {
    // Legion can interop with MPI in blocking and non-blocking
    // ways. You can use the calls to 'legion_wait_on_mpi' and
    // 'legion_handoff_to_mpi' in the same way as the MPI thread
    // does. Alternatively, you can get a phase barrier associated
    // with a LegionMPIHandshake object which will allow you to
    // continue launching more sub-tasks without blocking. 
    // For deferred execution we prefer the later style, but
    // both will work correctly.
    if (i < (total_iterations/2))
    {
      // This is the blocking way of using handshakes, it
      // is not the ideal way, but it works correctly
      // Wait for MPI to give us control to run our worker
      // This is a blocking call
      handshake.legion_wait_on_mpi();
      // Launch our worker task
      TaskLauncher worker_launcher(WORKER_TASK_ID, TaskArgument(NULL,0));
      Future f = runtime->execute_task(ctx, worker_launcher);
      // Have to wait for the result before signaling MPI
      f.get_void_result();
      // Perform a non-blocking call to signal
      // MPI that we are giving it control back
      handshake.legion_handoff_to_mpi();
    }
    else
    {
      // This is the preferred way of using handshakes in Legion
      TaskLauncher worker_launcher(WORKER_TASK_ID, TaskArgument(NULL,0));
      // We can user our handshake as a phase barrier
      // Record that we will wait on this handshake
      worker_launcher.add_wait_handshake(handshake);
      // Advance the handshake to the next version
      handshake.advance_legion_handshake();
      // Then record that we will arrive on this versions
      worker_launcher.add_arrival_handshake(handshake);
      // Launch our worker task
      // No need to wait for anything
      runtime->execute_task(ctx, worker_launcher);
    }
  }
}

void top_level_task(const Task *task, 
                    const std::vector<PhysicalRegion> &regions,
                    Context ctx, Runtime *runtime)
{
  printf("Hello from Legion Top-Level Task\n");
  // Both the application and Legion mappers have access to
  // the mappings between MPI Ranks and Legion address spaces
  // The reverse mapping goes the other way
  const std::map<int,AddressSpace> &forward_mapping = 
    runtime->find_forward_MPI_mapping();
  for (std::map<int,AddressSpace>::const_iterator it = 
        forward_mapping.begin(); it != forward_mapping.end(); it++)
    printf("MPI Rank %d maps to Legion Address Space %d\n", 
            it->first, it->second);

  int rank = -1, size = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Do a must epoch launch to align with the number of MPI ranks
  MustEpochLauncher must_epoch_launcher;
  LegionRuntime::Arrays::Rect<1> launch_bounds( LegionRuntime::Arrays::Point<1>(0), LegionRuntime::Arrays::Point<1>(size - 1));
  Domain launch_domain = Domain::from_rect<1>(launch_bounds);
  ArgumentMap args_map;
  IndexLauncher index_launcher(MPI_INTEROP_TASK_ID, launch_domain, 
                               TaskArgument(NULL, 0), args_map);
  must_epoch_launcher.add_index_task(index_launcher);
  runtime->execute_must_epoch(ctx, must_epoch_launcher);
}



int main(int argc, char *argv[])
{
   int numRanks = 1;
   int myRank = 0;

#if defined(COEVP_MPI)
#ifdef GASNET_CONDUIT_MPI
  // The GASNet MPI conduit requires special start-up
  // in order to handle MPI calls from multiple threads
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  // If you fail this assertion, then your version of MPI
  // does not support calls from multiple threads and you 
  // cannot use the GASNet MPI conduit
  if (provided < MPI_THREAD_MULTIPLE)
    printf("ERROR: Your implementation of MPI does not support "
           "MPI_THREAD_MULTIPLE which is required for use of the "
           "GASNet MPI conduit with the Legion-MPI Interop!\n");
  assert(provided == MPI_THREAD_MULTIPLE);
#else
  // Perform MPI start-up like normal for most GASNet conduits
  MPI_Init(&argc, &argv);
#endif

   MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
#endif

   
   
  int  timer = 0;
  int  sampling = 0;              //  By default, use adaptive sampling (but compiled in)
  int  redising = 0;              //  By default, do not use REDIS for database
  int  posixing = 0;              // POSIX Key/Value store
  int  hioing = 0;                //  HIO global store
  int  global_ns = 0;             //  By default, do not use a global nearest neighbor
  int  flanning = 0;              //  By default, do not use FLANN for nearest neighbor search
  int  logging = 0;               //  By default, do not use FLANN for nearest neighbor search
  int  flann_n_trees = 1;         // Default can be overridden using command line
  int  flann_n_checks = 20;       // Default can be overridden using command line
  int  file_parts = 0;
  int  debug_topology = 0;
  int  visit_data_interval = 0;  // Set this to 0 to disable VisIt data writing
  int  distributed_redis = 0;
  char logdb[1024] = {0};        // host and port of logging databse (e.g. cn1:6379)
  int use_vpsc = 0;              // toggle VPSC as the fine-scale-Model
  double c_scaling = 1.0;
  int heightElems = 26;
  int edgeElems = 16;
  double domStopTime = 1.e-1;
  int simStopCycle = 0;
  
  Lulesh luleshSystem;

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
  addArg("use_vpsc", 'm', 0, 'i',  &(use_vpsc),            0, "use VPSC fine scale model");
  addArg("convergence scaling", 'C', 1, 'd',  &(c_scaling), 0, "scaling factor for VPSC convergence"); 
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


  // Configure the Legion runtime with the rank of this process
  Runtime::configure_MPI_interoperability(myRank);
  // Register our task variants
  {
    TaskVariantRegistrar top_level_registrar(TOP_LEVEL_TASK_ID);
    top_level_registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    Runtime::preregister_task_variant<top_level_task>(top_level_registrar, 
                                                      "Top Level Task");
    Runtime::set_top_level_task_id(TOP_LEVEL_TASK_ID);
  }
  {
    TaskVariantRegistrar mpi_interop_registrar(MPI_INTEROP_TASK_ID);
    mpi_interop_registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    Runtime::preregister_task_variant<mpi_interop_task>(mpi_interop_registrar,
                                                        "MPI Interop Task");
  }
  {
    TaskVariantRegistrar worker_task_registrar(WORKER_TASK_ID);
    worker_task_registrar.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    Runtime::preregister_task_variant<worker_task>(worker_task_registrar,
                                                   "Worker Task");
  }
  // Create a handshake for passing control between Legion and MPI
  // Indicate that MPI has initial control and that there is one
  // participant on each side
  handshake = Runtime::create_handshake(true/*MPI initial control*/,
                                        1/*MPI participants*/,
                                        1/*Legion participants*/);
  // Start the Legion runtime in background mode
  // This call will return immediately
  Runtime::start(argc, argv, true/*background*/);
  // Run your MPI program like normal
  // If you want strict bulk-synchronous execution include
  // the barriers protected by this variable, otherwise
  // you can elide them, they are not required for correctness
  const bool strict_bulk_synchronous_execution = true;
  
  for (int i = 0; i < total_iterations; i++)
  {
    printf("MPI Doing Work on rank %d\n", myRank);
    if (strict_bulk_synchronous_execution)
      MPI_Barrier(MPI_COMM_WORLD);
    // Perform a handoff to Legion, this call is
    // asynchronous and will return immediately
    handshake.mpi_handoff_to_legion();
    // You can put additional work in here if you like
    // but it may interfere with Legion work

    // Wait for Legion to hand control back,
    // This call will block until a Legion task
    // running in this same process hands control back
    handshake.mpi_wait_on_legion();
    if (strict_bulk_synchronous_execution)
      MPI_Barrier(MPI_COMM_WORLD);
  }
  // When you're done wait for the Legion runtime to shutdown
  Runtime::wait_for_shutdown();
#ifndef GASNET_CONDUIT_MPI
  // Then finalize MPI like normal
  // Exception for the MPI conduit which does its own finalization
  MPI_Finalize();
#endif
  return 0;
  

  
  
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
  if (use_vpsc) {
   printf("Using VPSC fine-scale model\n");
  } else {
   printf("Using Taylor fine-scale model\n");
  }
  freeArgs();

/*
  if(timer)
  {
	// open output filestream here to avoid overhead
	// only myRank == 0 will do timing
	if(myRank==0)
	{
		//luleshSystem.timer = timer;
		luleshSystem.timerfile.open("timer.file");
		if(luleshSystem.timerfile.is_open())
		{
			for(int i=0;i<argc;i++)
			{
				luleshSystem.timerfile << argv[i] << " ";
			}
			
			luleshSystem.timerfile << std::endl;
			
	  	}
		else
		{
			std::cout << "Could not open timer.file" << std::endl;
		}
	}

  }
*/

   
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
  luleshSystem.ConstructFineScaleModel(sampling,global_modelDB,global_ann,flanning,flann_n_trees,flann_n_checks,global_ns,use_vpsc, c_scaling);
  
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
