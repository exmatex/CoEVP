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
   
   size_t sz = buf.size();
   
   std::string skey=uint128_to_string(key);
   hrc = hio_element_open(dataset, &element, skey,  HIO_FLAG_CREAT|HIO_FLAG_WRITE|HIO_FLAG_READ);
   hio_return(hrc);
   
   hrc = hio_element_write(element, 0, 0, &buf, sz, sz*sizeof(buf[0]));
   hio_return(hrc);
   
   //printf("Force element flush\n");
   //hrc = hio_element_flush(element, HIO_FLUSH_MODE_COMPLETE);
   //hio_return(hrc);

   hrc = hio_element_close(&element);
   hio_return(hrc);   
   
}

void  SingletonDB_HIO::erase(const uint128_t &key){
 
}

std::vector<double> SingletonDB_HIO::pull(const uint128_t &key) {
   printf("Read back data from element\n");
   
   std::string skey=uint128_to_string(key);
   size_t sz;
  
   hrc = hio_element_open(dataset, &element, skey,  HIO_FLAG_CREAT|HIO_FLAG_WRITE|HIO_FLAG_READ);
   hio_return(hrc);

   hrc = hio_element_size(element, &sz);
   std::vector<double> packedContainer(sz);
   hrc = hio_element_read(element, 0, 0, &packedContainer, sz, sz*sizeof(packedContainer[0]));
   hio_return(hrc);
   
   hrc = hio_element_close(&element);
   hio_return(hrc);
   

   return packedContainer;
}

std::vector<double> SingletonDB_HIO::pull_key(const uint128_t &key) {
   std::cout << "HIO DB: pull_key not implemented, contact vernon@lanl.gov" << std::endl;
   exit();
}


void  SingletonDB_HIO::check_hio_return(hio_return_t hrc)
{
   if(hrc==HIO_SUCCESS)
   {
       printf("HIO action succeeded: %d\n", hrc);
   }
   else
   {
       printf("HIO action failed: %d\n", hrc);
   }
}

SingletonDB_HIO::SingletonDB_HIO() {
   std::cout << "Create active HIO context..." << std::endl;
   hrc = hio_init_single(&hio_context, HIO_CONFIG, "", HIO_NAMESPACE);
   hio_return(hrc);
   
   std::cout << "Create/Open dataset instance" << std::endl;
   hrc  = hio_dataset_open(hio_context, &dataset, "kriging_database", 0, HIO_FLAG_CREAT|HIO_FLAG_WRITE|HIO_FLAG_READ, 1);
   hio_return(hrc);

}


SingletonDB_HIO::~SingletonDB_HIO() {
  std::cout << "Finishing HIO in SingletonDB_HIO..." << std::endl;
  
   printf("Close dataset \n");
   hrc = hio_dataset_close(&dataset);
   hio_return(hrc);
   printf("Destrony HIO context \n");   
   hio_fini(&hio_context);

}

#endif
