#include "domain.h"
#include "ElastoViscoPlasticity.h"

ConstitutiveData  shim_advance(Domain &domain, int k) {
  ConstitutiveData cm_data = domain.cm(k)->advance(domain.deltatime(),
                                                   domain.cm_vel_grad(k),
                                                   domain.cm_vol_chng(k),
                                                   domain.cm_state(k));
  return cm_data;
}
