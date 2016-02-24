#ifndef included_SingletonDB_HIO_h
#define included_SingletonDB_HIO_h

#include "SingletonDB_Backend.h"

#include <vector>
#include <hio.h>
#define uint128_t unsigned __int128

class SingletonDB_HIO : public SingletonDB_Backend{
 public:
 
  SingletonDB_HIO();
  ~SingletonDB_HIO();

  virtual void  push(const uint128_t &key, const std::vector<double>& buf, const unsigned long key_length);
  virtual void  erase(const uint128_t &key);
  virtual std::vector<double> pull(const uint128_t &key);
  virtual std::vector<double> pull_key(const uint128_t &key);

private:
  int hrc;
  hio_context_t hio;
  std::string prefix;
  enum dataset_state { HREAD, HWRITE };
  hio_dataset_t dataset;
  dataset_state dstate;

  void  check_hio_return(int hrc, hio_context_t hio, char description[])
  {
     if(hrc==HIO_SUCCESS)
     {
		(void) hio_err_print_all(hio, stderr, description);
     }
//     else
//     {
//         printf("HIO action failed: %d\n", hrc);
//     }

  }

};


#endif // included_SingletonDB_HIO_h
