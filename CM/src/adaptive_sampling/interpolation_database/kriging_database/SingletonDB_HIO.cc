#ifdef HIO

#include "SingletonDB_HIO.h"

#ifndef HIO_NAMESPACE
#define HIO_NAMESPACE "CoEVP"
#endif
#ifndef HIO_CONFIG
#define HIO_CONFIG "config.hio"
#endif

#include <iostream>
#include <string>
#include <assert.h>
#include <stdexcept>
#include <unistd.h>
#include <hio.h>


static std::string uint128_to_string(const uint128_t &in){
   uint64_t *in64 = (uint64_t *)&in; 
   return std::to_string(*in64)+std::to_string(*(in64+1));
}
      
//  Will eventually be something like add_points
void  SingletonDB_HIO::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {

   hio_element_t element;
   
//   size_t sz;
   int64_t sz = buf.size();
   printf("IN: Number of elements %d\n", buf.size());
   printf("Element[0] %g\n", buf[0]);
//   for (int i = sz - 1; i >= 0; i--) 
//     std::cout << buf[0];

   printf("Create New HIO Element\n");
   std::string skey=uint128_to_string(key);
   hrc = hio_element_open(dataset, &element, skey.c_str(),  HIO_FLAG_CREAT|HIO_FLAG_WRITE);
   check_hio_return(hrc);
   
   printf("Write Key to Element\n");
   hrc = hio_element_write(element, 0, 0, &buf, sz, sizeof(double));
   check_hio_return(hrc-sz*sizeof(double));
   
   printf("Force element flush\n");
//   hrc = hio_element_flush(element, HIO_FLUSH_MODE_LOCAL);
   hrc = hio_dataset_flush(dataset, HIO_FLUSH_MODE_LOCAL);
   check_hio_return(hrc);

   printf("Close Element\n");
   hrc = hio_element_close(&element);
   check_hio_return(hrc);   


   
}

void  SingletonDB_HIO::erase(const uint128_t &key){
   std::cout << "HIO DB: erase not implemented, contact vernon@lanl.gov" << std::endl;
   exit(1);
 
}

std::vector<double> SingletonDB_HIO::pull(const uint128_t &key) {


   hio_element_t element;
   
   std::string skey=uint128_to_string(key);
//   size_t sz;
   int64_t sz;

   printf("Open Element for Reading\n");
  
   hrc = hio_element_open(dataset, &element, skey.c_str(),  HIO_FLAG_READ);
   check_hio_return(hrc);

   printf("Get Size of Element\n");
   hrc = hio_e_size(element, &sz);
   check_hio_return(hrc);
   sz /= sizeof(double);
   printf("OUT: Size of element %d\n", sz);
   std::vector<double> packedContainer(sz);

   printf("Read Element Contents\n");
   hrc = hio_element_read(element, 0, 0, &packedContainer, sz, sizeof(double));
   check_hio_return(hrc);
   printf("Element[0] %g\n", packedContainer[0]);
//   for (int i = sz - 1; i >= 0; i--) 
//     std::cout << packedContainer[i];

   printf("Size of packedContainer %d\n", packedContainer.size());

   
   printf("Close Element\n");

   hrc = hio_element_close(&element);
   check_hio_return(hrc);
   
   return packedContainer;
}

std::vector<double> SingletonDB_HIO::pull_key(const uint128_t &key) {
   std::cout << "HIO DB: pull_key not implemented, contact vernon@lanl.gov" << std::endl;
   exit(1);
}



SingletonDB_HIO::SingletonDB_HIO() {
   std::cout << "HIO API VERSION: " << HIO_API_VERSION << std::endl;
   std::cout << "Create active HIO context."  << std::endl;
   hrc = hio_init_single(&hio, HIO_CONFIG, "", HIO_NAMESPACE);
   check_hio_return(hrc);
   
   std::cout << "Create/Open dataset instance" << std::endl;
   hrc  = hio_dataset_open(hio, &dataset, "kriging_database", 0, HIO_FLAG_CREAT|HIO_FLAG_WRITE|HIO_FLAG_READ, HIO_SET_ELEMENT_UNIQUE);
   check_hio_return(hrc);

}


SingletonDB_HIO::~SingletonDB_HIO() {
  std::cout << "Finishing HIO in SingletonDB_HIO..." << std::endl;
  
   printf("Close Dataset \n");
   hrc = hio_dataset_close(&dataset);
   check_hio_return(hrc);
   printf("Destroy HIO Context \n");   
   hio_fini(&hio);

}

#endif

