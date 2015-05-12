//
// File:        DBKrigingModelObject.h
// Package:     kriging coupler
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: Specialization of the DBModelObject to kriging.
//

#ifndef included_krigcpl_DBKrigingModelObject_h
#define included_krigcpl_DBKrigingModelObject_h

#include "DBModelObject.h"

#include <base/InterpolationModel.h>

namespace krigcpl {

    typedef DBModelObject<krigalg::InterpolationModelPtr> DBKrigingModelObject;

}
  


#endif // included_krigcpl_DBKrigingModelObject_h
