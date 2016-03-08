//  Small test code to experiment with libcircle.
//
//  mpic++ -std=c++11 -g -O3 lctest.cc -o lctest -I ../serverize/circle/include -I ../redis/hiredis ../redis/hiredis/libhiredis.a ../serverize/circle/lib/libcircle.a
//
//  mpirun -np 2 ./lctest  (remember to start a REDIS server on target node)
//
//  mpirun --map-by node -np 2 ./lctest   (fails...hangs?)

#include <iostream>
#include <stdexcept>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <libcircle.h>
#include <mpi.h>
#include "hiredis.h"

int timestep = 0;


//  Get functions parameters (unpack protobuf-encoded string), call
//  function, and write results to database.
void processTasks(CIRCLE_handle *handle) {
  char taskString[CIRCLE_MAX_STRING_LEN];
  handle->dequeue(&taskString[0]);
  //  "unpack"
  int  ts;
  int  k;
  sscanf(taskString, "%d %d", &ts, &k);

  int numRanks;
  int myRank;
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
  char my_node[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(my_node, &name_len);

  //  Fake putting results into REDIS.
  redisContext  *redis = redisConnect("cn2", 6379);   // [HACK] start manually on MPI node
  if (redis != NULL && redis->err) {
    throw std::runtime_error("Error connecting to redis for circle");
  }
  char key[100];
  sprintf(key, "%d:%d", ts, k);
  redisReply *reply = (redisReply *)redisCommand(redis, "SET %s %s", key, "results");
  if (!reply) {
    std::cerr << "No connection to redis for libcircle...continuing" << std::endl;
  }
  freeReplyObject(reply);
  redisFree(redis);
  
  // std::cout << my_node << "/" << myRank << "/" << numRanks << " task string: " << taskString << std::endl;
}


//  Setup (via protobuf) parameters to function call.  Would have an
//  enqueue for each of the 'k' cells in the mesh.
void buildTasks(CIRCLE_handle *handle) {
  char taskString[CIRCLE_MAX_STRING_LEN];
  for(int i=0; i<5; i++) {
    //  "pack"
    sprintf(taskString, "%d %d", timestep, i);
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
    sleep(1);
    std:: cout << "------" << std::endl;
    //  Fake getting results from REDIS.
    redisContext  *redis = redisConnect("cn2", 6379);   // [HACK]
    if (redis != NULL && redis->err) {
      throw std::runtime_error("Error connecting to redis for circle");
    }
    redisReply *reply = (redisReply *)redisCommand(redis, "KEYS *");
    if (!reply) {
      std::cerr << "No connection to redis for libcircle...continuing" << std::endl;
    }
    if (reply->type != REDIS_REPLY_ARRAY) {
      throw std::runtime_error("Wrong redis return type for libcircle KEYS *");
    }
    if(reply->elements < 1) {
      throw std::runtime_error("Number of redis reply elements wrong on KEYS *, got: "
                               + std::to_string(reply->elements));
    }
    std::cout << "TS: " << i << " (" << reply->elements << ")" << std::endl;
    for (int j=0; j<reply->elements; j++) {
      redisReply  *k_reply = (redisReply *)reply->element[j];
      if (k_reply->type != REDIS_REPLY_STRING) {
        throw std::runtime_error("Wrong redis return type for libcircle getting a key");
      }
      std::cout << k_reply->str << std::endl;
    }
    freeReplyObject(reply);   //  free the KEYS * results (and sub-results)
    reply = (redisReply *)redisCommand(redis, "FLUSHDB");
    if (!reply) {
      std::cerr << "No connection to redis for libcircle...continuing" << std::endl;
    }
    freeReplyObject(reply);   //  free the FLUSHDB command results
    redisFree(redis);
  }

  CIRCLE_finalize();
  return 0;
}
