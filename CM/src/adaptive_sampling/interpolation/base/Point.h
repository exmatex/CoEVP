//
// File:        Point.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Abstract base class for point used in Kriging computations.
//
// $Id: Point.h,v 1.1 2005/08/23 21:12:40 knap2 Exp $
//
// $Log: Point.h,v $
// Revision 1.1  2005/08/23 21:12:40  knap2
// Initial source.
//
//

#if !defined(included_krigalg_Point)
#define included_krigalg_Point

#ifndef included_krigalg_Vector
#include "Vector.h"
#endif


namespace krigalg {

  class Point : public Vector {

  public:
    //
    // construction/destruction
    //
    Point();
    Point(int n);
    Point(int n, double init);
    Point(int n, const double * coordinates);
    virtual ~Point();

    //
    // output
    //

    virtual std::ostream & print(std::ostream & outputStream) const;

  };


  //
  // addition
  //

  Point & operator+=(Point       & point1,
		     const Point & point2); // invalid - not implemented 
  Point & operator+=(Point        & point,
		     const Vector & vector);

  Point operator+(const Point & point1, const Point & point2);//not implemented
  Point operator+(const Point & point, const Vector & vector);
  Point operator+(const Vector & vector, const Point & point);

  //
  // subtraction
  //

  Vector operator-(const Point & point1, const Point & point2);

  //
  // distance function
  //

  double distance(const Point & point1, const Point & point2);

  //
  // output
  //

  std::ostream & operator<<(std::ostream & outputStream,
			    const Point  & point);

}

#ifndef DEBUG_NO_INLINE
#include "Point.I"
#endif // DEBUG_NO_INLINE

#endif // included_krigalg_Point
