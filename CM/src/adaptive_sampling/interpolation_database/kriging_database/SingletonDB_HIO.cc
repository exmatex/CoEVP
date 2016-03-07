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

#include "KeyToString.h"

void  SingletonDB_HIO::push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length) {

//   hio_context_t hio;
//  hio_return_t hrc;
//   int  hrc;

//   std::cout << "Open HIO context" << std::endl;
//   hrc = hio_init_single(&hio, HIO_CONFIG, "", HIO_NAMESPACE);
//   check_hio_return(hrc);

	switch(dstate)
    {
		case HREAD:
			hrc = hio_dataset_close(&dataset);
		    check_hio_return(hrc, hio, "HIO: Close Dataset Read Handle"); 
    		hrc  = hio_dataset_open(hio, &dataset, prefix.c_str(), 0, HIO_FLAG_WRITE, HIO_SET_ELEMENT_UNIQUE);
			check_hio_return(hrc, hio, "HIO: Open Dataset Write Handle");
			dstate = HWRITE;
			break;
	}

    hio_element_t element;

//   size_t sz;
    int64_t sz = buf.size();

    std::string skey=uint128_to_string(key);

    hrc = hio_element_open(dataset, &element, skey.c_str(),  HIO_FLAG_CREAT|HIO_FLAG_WRITE);
    check_hio_return(hrc, hio, "HIO: Create new element");

    hrc = hio_element_write(element, 0, 0, buf.data(), sz, sizeof(double));
    check_hio_return(hrc-sz*sizeof(double), hio, "HIO: Write buffer to element");
   
//   printf("Force element flush\n");
//   hrc = hio_element_flush(element, HIO_FLUSH_MODE_LOCAL);
//   hrc = hio_dataset_flush(dataset_writer, HIO_FLUSH_MODE_COMPLETE);
//   check_hio_return(hrc);

    hrc = hio_element_close(&element);
    check_hio_return(hrc, hio, "HIO: Close Element");   

   
//   std::cout << "Close HIO context" << std::endl;
//   hrc = hio_fini(&hio);
//   check_hio_return(hrc);
   
}

void  SingletonDB_HIO::erase(const uint128_t &key){
    std::cout << "HIO: erase not implemented, contact vernon@lanl.gov" << std::endl;
    exit(1);
 
}

std::vector<double> SingletonDB_HIO::pull(const uint128_t &key) {

//   hio_context_t hio;
//   int  hrc;

    switch(dstate)
    {
        case HWRITE:
            hrc = hio_dataset_close(&dataset);
            check_hio_return(hrc, hio, "HIO: Close Dataset Write Handle"); 
            hrc  = hio_dataset_open(hio, &dataset, prefix.c_str(), 0, HIO_FLAG_READ, HIO_SET_ELEMENT_UNIQUE);
            check_hio_return(hrc, hio, "HIO: Open Dataset Read Handle");
            dstate = HREAD;
            break;
    }


//   std::cout << "Open HIO context" << std::endl;
//   hrc = hio_init_single(&hio, HIO_CONFIG, "", HIO_NAMESPACE);
//   check_hio_return(hrc);




    hio_element_t element;
   
    std::string skey=uint128_to_string(key);
//   size_t sz;
    int64_t sz;

  
    hrc = hio_element_open(dataset, &element, skey.c_str(),  HIO_FLAG_READ);
    check_hio_return(hrc, hio, "HIO: Open Element for reading");

    hrc = hio_element_size(element, &sz);
    check_hio_return(hrc, hio, "HIO: Query element size");
    sz /= sizeof(double);
    std::vector<double> packedContainer(sz);

    hrc = hio_element_read(element, 0, 0, packedContainer.data(), sz+10, sizeof(double));
    check_hio_return(hrc-sz*sizeof(double), hio, "HIO: Read element to buffer");


    hrc = hio_element_close(&element);
    check_hio_return(hrc, hio, "HIO: Close Element");

   
    return packedContainer;
}

std::vector<double> SingletonDB_HIO::pull_key(const uint128_t &key) {
    std::cout << "HIO: pull_key not implemented, contact vernon@lanl.gov" << std::endl;
    exit(1);
}



SingletonDB_HIO::SingletonDB_HIO() {
    std::cout << "HIO API VERSION: " << HIO_API_VERSION << std::endl;
    hrc = hio_init_single(&hio, HIO_CONFIG, "", HIO_NAMESPACE);
    check_hio_return(hrc, hio, "HIO: Open HIO Context");   
   
// shameful prefix generation to deal with current HIO thread safety issues

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);
    prefix = std::string(hostname);
    prefix += std::to_string(getpid());
    prefix += "kriging_database";

	// create or open the dataset for the first time

    hrc  = hio_dataset_open(hio, &dataset, prefix.c_str(), 0, HIO_FLAG_CREAT|HIO_FLAG_WRITE, HIO_SET_ELEMENT_UNIQUE);
    check_hio_return(hrc, hio, "HIO: Create Dataset");
    if(hrc != HIO_SUCCESS)
    {
    	hrc  = hio_dataset_open(hio, &dataset, prefix.c_str(), 0, HIO_FLAG_WRITE, HIO_SET_ELEMENT_UNIQUE);
		check_hio_return(hrc, hio, "HIO: Open existing Dataset");
		if(hrc != HIO_SUCCESS)
		{
			std::perror("HIO: Can't create or open existing dataset. Check HIO lib compatibility?\n");
			exit(1);

		}
	}
	dstate = HWRITE;
}


SingletonDB_HIO::~SingletonDB_HIO() {
    std::cout << "HIO: Finishing HIO in SingletonDB_HIO..." << std::endl;
    hio_fini(&hio);

}

#endif

