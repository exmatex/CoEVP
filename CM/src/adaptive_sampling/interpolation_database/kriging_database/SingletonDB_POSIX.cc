#include "SingletonDB_POSIX.h"

#ifndef POSIX_PATH
#define POSIX_PATH "."

#include <iostream>
#include <string>
#include <assert.h>
#include <stdexcept>
#include <unistd.h>
#include <fstream>
#include <cstdio>

static std::string uint128_to_string(const uint128_t &in){
   uint64_t *in64 = (uint64_t *)&in; 
   return std::to_string(*in64)+std::to_string(*(in64+1));
}
      

void  SingletonDB_POSIX::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {
  size_t sz = buf.size();
  std::string skey=uint128_to_string(key);
  strcat(skey, dbPath);
  std::ofstream keyfile(skey, std::ios::out | std::ios::binary);
  keyfile.write((char *) &buf[0], sz*sizeof(buf[0]));
  keyfile.close();
}

void  SingletonDB_POSIX::erase(const uint128_t &key){
  std::string skey=uint128_to_string(key);
  strcat(skey, dbPath);
  int reply = std::remove(skey.c_str());
}

std::vector<double> SingletonDB_POSIX::pull(const uint128_t &key) {
  std::string skey=uint128_to_string(key);
  strcat(skey, dbPath);
  std::ifstream keyfile(skey, std::ios::in | std::ios::binary);
  keyfile.seekg(0, keyfile.end);
  size_t sz = keyfile.tellg();
  keyfile.seekg(0, keyfile.beg);
  std::vector<double> packedContainer(sz);
  keyfile.read((char *) &packedContainer[0], sz*sizeof(double));

  return packedContainer;
}

std::vector<double> SingletonDB_POSIX::pull_key(const uint128_t &key) {
   std::cout << "POSIX DB: pull_key not implemented, contact vernon@lanl.gov" << std::endl;
   exit();
  //std::string skey=uint128_to_string(key);
  //strcat(skey, dbPath);
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
  std::cout << "POSIX DB, no initialization needed" << std::endl;
  dbPath = POSIX_PATH;
}

//  Shutdown redis databse and print some simple info about the accumulated database.
//
SingletonDB_POSIX::~SingletonDB_POSIX() {
  std::cout << "POSIX DB, nothing to destroy..." << std::endl;
}

