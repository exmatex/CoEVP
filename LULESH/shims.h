#ifndef SHIMS_H_
#define SHIMS_H_

class  Domain;
class  ConstitutiveData;

struct WrapReturn {
  ConstitutiveData  *cm_data;
  void              *state;
};
  
struct WrapReturn  *wrap_advance(Domain &, int);

#endif

