#ifndef PACK_H_
#define PACK_H_

class  Domain;
class  Tensor2Gen;
class  ConstitutiveData;

std::string  unpackAdvanceAndCall(Domain &, std::string, int *, ConstitutiveData *, void *);
void  unpackAdvance(std::string, int *, ConstitutiveData *, void *);
std::string  packAdvance(int, ConstitutiveData &, double, Tensor2Gen &, double, void *);
#endif

