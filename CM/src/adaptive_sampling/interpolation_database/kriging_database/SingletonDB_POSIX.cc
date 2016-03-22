#include "SingletonDB_POSIX.h"
#include <sys/stat.h>

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <assert.h>
#include <stdexcept>
#include <unistd.h>
#include <fstream>
#include <cstdio>
//#include <unistd.h>

#include "KeyToString.h"


inline bool file_exists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}


void  SingletonDB_POSIX::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {
  unsigned long sz = buf.size();
  // generate path to file
  std::string skey=uint128_to_string(key);
  std::string skey_path = dbPath + "/" + skey;
  // check key does not exist
  std::ifstream infile(skey_path);
  if(infile.good())
  { 
    //std::cout << "Key Hit, returning" << std::endl;
    infile.close();
	return;
  } 

  // write with prefix
//  std::cout << "Key Miss, writing" << std::endl;

  std::string prefix_skey = dbPath + "/" + prefix + skey;  
  std::ofstream keyfile(prefix_skey, std::ios::out | std::ios::binary);

//  keyfile.write((char *) sz, sizeof(unsigned long));
//  keyfile.write((char *) key_length, sizeof(unsigned long));
  keyfile.write((char *) buf.data(), sz*sizeof(buf[0]));
  keyfile.close();


// std::cout << "Key written, renaming" << std::endl;
 
//  printf("Size of buf %d\n", buf.size() * sizeof(buf[0]));

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

  std::vector<double> packedContainer(sz+2);
//  keyfile.read((char *) &packedContainer[0], sz*sizeof(double));
  keyfile.read((char *) packedContainer.data(), sz);

//  exit(1);
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
  std::ifstream configFile;
  configFile.open("config.posix", std::ifstream::in);
  if(configFile.good() == true)
  {
      ///TODO: Make this remotely safe. Possibly by using the boost json parser we already need anyway
      std::string configLine;
      getline(configFile, configLine);
      //Remove "data_roots ="
      configLine.erase(0, 12);
      //If there is a leading space, get rid of it
      if(configLine[0] == ' ')
      {
       configLine.erase(0, 1);
      }
      //And now just copy the string
      dbPath = configLine;
      configFile.close();
  }
  else
  {
      dbPath = "posix";
  }

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

