//
// File:        Vector.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing MTL-aware vector.
//
// $Id$
//
// $Log$
//

#if !defined(included_krigalg_Vector)
#define included_krigalg_Vector

#ifndef included_iosfwd
#define included_iosfwd
#include <iosfwd>
using namespace std;
#endif

#ifndef included_mtl_dense1D
#define included_mtl_dense1D
#include <mtl/dense1D.h>
#endif

namespace krigalg {

  // typedef mtl::dense1D<double> Vector;

  class Vector : public mtl::dense1D<double> {

  public:
    //
    // construction/destruction
    //
    Vector();
    Vector(int n);
    Vector(int n, double init);
    Vector(int n, const double * coordinates);
    virtual ~Vector();

    //
    // output
    //

    virtual std::ostream & print(std::ostream & outputStream) const;

  };

  //
  // addition
  //

  Vector & operator+=(Vector       & vector1,
		      const Vector & vector2);
  Vector operator+(const Vector & vector1,
		   const Vector & vector2);

  //
  // subtraction
  //

  Vector & operator-=(Vector       & vector1,
		      const Vector & vector2);
  Vector operator-(const Vector & vector1,
		   const Vector & vector2);

  //
  // scalar product
  //

  Vector & operator*=(Vector & vector1,
		      double alpha);
  Vector operator*(const Vector & vector1,
		   double alpha);


  //
  // dot product
  //

  double dot(const Vector & vector1,
	     const Vector & vector2);

  //
  // output
  //

  std::ostream & operator<<(std::ostream & outputStream,
			    const Vector & vector);
  
}

#ifndef DEBUG_NO_INLINE
#include "Vector.I"
#endif // DEBUG_NO_INLINE

#endif // included_krigalg_Vector
