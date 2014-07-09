//
// File:        SecondMoment.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Implement operations related to computing tensor of inertia
// Description: for a collection of Points. The tensor is computed wrt the 
// Description: center of mass.
//

#if !defined(included_krigalg_SecondMoment)
#define included_krigalg_SecondMoment

#ifndef included_krigalg_Point
#include <base/Point.h>
#endif

#ifndef included_krigalg_Matrix
#include <base/Matrix.h>
#endif

#ifndef included_krigalg_Vector
#include <base/Vector.h>
#endif

#ifndef included_iosfwd
#define included_iosfwd
#include <iosfwd>
using namespace std;
#endif

#ifndef included_vector
#define included_vector
#include <vector>
using namespace std;
#endif

#ifndef included_utility
#define included_utility
#include <utility>
using namespace std;
#endif

namespace krigalg {

  //
  // forward declarations
  //

  class SecondMoment;

  //
  //
  //

  class SecondMoment {

  public:
    //
    // construction/destruction
    //
    SecondMoment(const std::vector<Point> & points);
    ~SecondMoment();

    //
    // get handle to center of mass
    //

    const Point & getCenterMass() const;

    //
    // compute eigensystem
    //

    std::pair<Vector, Matrix> computeEigenSystem() const;

    //
    // output
    //
    
    friend std::ostream & operator<<(std::ostream       & outputStream,
				     const SecondMoment & inertiaTensor);
    
  private:
    //
    // copy construction/assignment
    //
    SecondMoment(const SecondMoment &);
    const SecondMoment & operator=(const SecondMoment &);

    //
    // data
    //

  public:

  private:
    Point   _centerMass;
    Matrix  _data;       // FIXME: needs to be using symmetric storage

  };

}

#ifndef DEBUG_NO_INLINE
#include "SecondMoment.I"
#endif // DEBUG_NO_INLINE

#endif // included_krigalg_SecondMoment
