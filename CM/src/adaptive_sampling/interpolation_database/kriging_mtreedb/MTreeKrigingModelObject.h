//
// File:        MTreeKrigingModelObject.h
// Package:     kriging coupler
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: Specialization of the MTreeModelObject to kriging.
//

#ifndef included_krigcpl_MTreeKrigingModelObject_h
#define included_krigcpl_MTreeKrigingModelObject_h

#include "base/MTreeModelObject.h"

#include <base/InterpolationModel.h>

namespace krigcpl {

    typedef MTreeModelObject<krigalg::InterpolationModelPtr> MTreeKrigingModelObject;

}
  


#endif // included_krigcpl_MTreeKrigingModelObject_h
