//
// File:        Value.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Abstract base class for value used in Kriging computations.
//
// $Id: Value.h,v 1.1 2005/08/23 21:12:40 knap2 Exp $
//
// $Log: Value.h,v $
// Revision 1.1  2005/08/23 21:12:40  knap2
// Initial source.
//
//

#if !defined(included_krigalg_Value)
#define included_krigalg_Value

#ifndef included_krigalg_Vector
#include "Vector.h"
#endif


namespace krigalg {

  class Value : public Vector {

  public:
    //
    // construction/destruction
    //
    Value();
    Value(int n);
    Value(int n, double init);
    Value(int n, const double * components);
    Value(const Vector & vector);
    virtual ~Value();

    //
    // output
    //

    virtual std::ostream & print(std::ostream & outputStream) const;
    
  private:

  };

  //
  // output
  //

  std::ostream & operator<<(std::ostream & outputStream,
			    const Value  & value);

}

#endif // included_krigalg_Value
