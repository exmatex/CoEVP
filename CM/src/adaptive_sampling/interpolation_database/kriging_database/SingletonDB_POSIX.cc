#include "SingletonDB_POSIX.h"
#include <sys/stat.h>

#ifndef POSIX_PATH
#define POSIX_PATH "posix"
#endif

#include <iostream>
#include <string>
#include <assert.h>
#include <stdexcept>
#include <unistd.h>
#include <fstream>
#include <cstdio>
//#include <unistd.h>

static std::string uint128_to_string(const uint128_t &in){
   uint64_t *in64 = (uint64_t *)&in; 
   return std::to_string(*in64)+std::to_string(*(in64+1));
}
      

void  SingletonDB_POSIX::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {
  size_t sz = buf.size();
  // generate path to file
  std::string skey=uint128_to_string(key);
  std::string skey_path = dbPath + "/" + skey;
  // check key does not exist
  std::ifstream infile(skey_path);
  if(infile.good())
  { 
//    std::cout << "Key Hit, returning" << std::endl;
    infile.close();
	return;
  }
  // write with prefix
//  std::cout << "Key Miss, writing" << std::endl;

  std::string prefix_skey = dbPath + "/" + prefix + skey;  
  std::ofstream keyfile(prefix_skey, std::ios::out | std::ios::binary);
  keyfile.write((char *) &buf[0], sz*sizeof(buf[0]));
  keyfile.close();


// std::cout << "Key written, renaming" << std::endl;
 
  int rc = std::rename(prefix_skey.c_str(), skey_path.c_str()); 
  if(rc) { std::perror("Key Write Soft Collision");}
}

void  SingletonDB_POSIX::erase(const uint128_t &key){
//  std::string skey=uint128_to_string(key);
//  skey = dbPath + "/" + skey;
//  int reply = std::remove(skey.c_str());
   std::cout << "POSIX DB: erase not officially implemented, contact vernon@lanl.gov" << std::endl;
   exit(1);

}

std::vector<double> SingletonDB_POSIX::pull(const uint128_t &key) {
  std::string skey=uint128_to_string(key);
  skey = dbPath + "/" + skey;
  std::ifstream keyfile(skey, std::ios::in | std::ios::binary);
  keyfile.seekg(0, keyfile.end);
  size_t sz = keyfile.tellg();
  keyfile.seekg(0, keyfile.beg);
  std::vector<double> packedContainer(sz);
  keyfile.read((char *) &packedContainer[0], sz*sizeof(double));

  return packedContainer;
}

std::vector<double> SingletonDB_POSIX::pull_key(const uint128_t &key) {
   std::cout << "POSIX DB: pull_key not officially implemented, contact vernon@lanl.gov" << std::endl;
   exit(1);
  //std::string skey=uint128_to_string(key);
  //skey = dbPath + "/" + skey;
  //std::ifstream keyfile(skey, std::ios::in | std::ios::binary);
  //keyfile.seekg(0, keyfile.end);
  //size_t sz = keyfile.tellg();
  //keyfile.seekg(0, keyfile.beg);
  //
  //std::vector<double> packedContainer(sz);
  //keyfile.read((char *) &packedContainer[0], sz*sizeof(double));
  //
  //return packedContainer;
}


SingletonDB_POSIX::SingletonDB_POSIX() {
  std::cout << "POSIX DB, validating/creating path:" << std::endl;
  dbPath = POSIX_PATH;

  struct stat st = {0};
  if (stat(dbPath.c_str(), &st) == -1) {
      std::cout << "POSIX DB, path does not exist - creating..." << std::endl;
      mkdir(dbPath.c_str(), 0700);
  }
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  prefix = std::string(hostname);
  prefix += std::to_string(getpid());

}

//  Shutdown redis databse and print some simple info about the accumulated database.
//
SingletonDB_POSIX::~SingletonDB_POSIX() {
  std::cout << "POSIX DB, nothing to destroy..." << std::endl;
}

