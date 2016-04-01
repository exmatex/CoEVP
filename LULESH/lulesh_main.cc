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

#if defined(PROTOBUF)
#include <libcircle.h>
#include <mpi.h>
#include <unistd.h>
#include "hiredis.h"
//  Hacky globals necessary for static libcircle callbacks to get handle to advance() function.
Domain *circleDomain;
redisContext *circleDB;
#endif

//  Command line option parsing (using Sriram code from old days)
#include "cmdLineParser.h"

int main(int argc, char *argv[])
{
   int numRanks = 1;
   int myRank = 0;
#if defined(COEVP_MPI)
   MPI_Init(&argc, &argv) ;
   MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
#endif
   //  Can't do MPI and libcircle at the same time--libcircle deals with
   //  MPI in it's own way.
#if defined(PROTOBUF)
   //  We're NOT doing domain decomposition with MPI--just libcircle.
   //  so, as far as lulesh knows, there is only one rank. The libcircle
   //  code exclusively uses the other ranks--no need for lulesh to even
   //  know about that.
   int bogusRank = CIRCLE_init(argc, argv, CIRCLE_DEFAULT_FLAGS);
   CIRCLE_cb_process(Lulesh::processCircleTasks);
   CIRCLE_enable_logging(CIRCLE_LOG_ERR);
   MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   std::cout << "Initialized libcircle..." << std::endl;
#endif
  
  int  sampling = 0;              //  By default, use adaptive sampling (but compiled in)
  int  redising = 0;              //  By default, do not use REDIS for database
  int  posixing = 0;              // POSIX Key/Value store
  int  hioing = 0;                //  HIO global store
  int  global_ns = 0;             //  By default, do not use a global earest neighbor
  int  flanning = 0;              //  By default, do not use FLANN for nearest neighbor search
  int  logging = 0;               //  By default, don't do logging
  int  circling = 0;              //  By default, don't use libcircle
  int  flann_n_trees = 1;         // Default can be overridden using command line
  int  flann_n_checks = 20;       // Default can be overridden using command line
  int  file_parts = 0;
  int  debug_topology = 0;
  int  visit_data_interval = 0;  // Set this to 0 to disable VisIt data writing
  int  distributed_redis = 0;
  char logdb[1024] = {0};         // host and port of logging databse (e.g. cn1:6379)
  char circledb[1024] = {0};       // host and port of libcircle results databse (e.g. cn1:6379)
  int heightElems = 26;
  int edgeElems = 16;
  double domStopTime = 1.e-1;
  
  Lulesh luleshSystem;

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
  addArg("Circle",      'C', 1, 's',  &(circledb),   sizeof(circledb), "libcircle results to REDIS at hostname:port");
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
  if (strlen(circledb) != 0) {
    std::cout << "circling to: " << circledb << std::endl;
    circling = 1;
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
#if defined(COEVP_MPI) || defined(PROTOBUF)
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

#if defined(PROTOBUF)
   if (circling) {
     circleDomain = &luleshSystem.domain;  // Set global pointer to this instance's Domain
     //  Move all this out to resource manager at some point.
     //  Separate the hostname from the port (input must be of form "host:nnn")
     std::vector<std::string> elems;
     std::stringstream ss(circledb);
     std::string item;
     while (std::getline(ss, item, ':')) {
       elems.push_back(item);
     }
     if (elems.size() != 2) {
       throw std::runtime_error("Invalid circledb host:port (" + std::string(circledb) + ") specified");
     }
     //  Try connecting to the libcircle database
     circleDB = redisConnect(elems[0].c_str(), std::stoi(elems[1]));
     if (circleDB != NULL && circleDB->err) {
       throw std::runtime_error("Error connecting to redis for libcircle, please start one on "
                                + std::string(circledb));
     }
     int circleRank;
     MPI_Comm_rank(MPI_COMM_WORLD, &circleRank);
     if (circleRank == 0) {
       redisReply *reply = (redisReply *) redisCommand(circleDB, "DBSIZE");
       if (!reply) {
         throw std::runtime_error("No connection to redis server for libcircle, please start one on "
                                  + std::string(circledb));
       }
       if (reply->type != REDIS_REPLY_INTEGER){
         throw std::runtime_error("Wrong redis return type (when opening for libcircle)");
       }
       if(reply->integer > 0) {
         std::cout << "...existing database found (" << reply->integer << " keys), erasing it" << std::endl;
         redisCommand(circleDB, "FLUSHDB");    // Never supposed to fail :-)
       }
       freeReplyObject(reply);
     }
     std::cout << "...connected for libcircle" << std::endl;
   }
#endif
   
  // Construct fine scale models
  luleshSystem.ConstructFineScaleModel(sampling,global_modelDB,global_ann,flanning,flann_n_trees,flann_n_checks,global_ns);
  
  // Exchange nodal mass
  luleshSystem.ExchangeNodalMass();

  // Simulate 
#if defined(PROTOBUF)
   if (circling) {
     luleshSystem.gogo(myRank,numRanks,sampling,visit_data_interval,file_parts,debug_topology);
   } else {
     luleshSystem.go(myRank,numRanks,sampling,visit_data_interval,file_parts,debug_topology);
   }
#endif
   
  // Only do this is we have actually opened a REDIS connection.
#if defined(LOGGER)
  if (logging) {
    delete(&logger); 
  }
#endif
  
#if defined(PROTOBUF)
  if (circling) {
    redisFree(circleDB);
  }
  CIRCLE_finalize();
#endif

#if defined(COEVP_MPI)
   MPI_Finalize() ;
#endif

  return 0;
}
