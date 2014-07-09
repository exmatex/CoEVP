//
// File:        Dimension.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Base class for the concept of dimension.
//
// $Id$
//
// $Log$
//

#if !defined(included_krigalg_Dimension)
#define included_krigalg_Dimension

#ifndef included_mtl_dimension
#define included_mtl_dimension
#include <mtl/dimension.h>
#endif

namespace krigalg {

  typedef mtl::dimension<int> Dimension;

}

#endif // included_krigalg_Dimension
