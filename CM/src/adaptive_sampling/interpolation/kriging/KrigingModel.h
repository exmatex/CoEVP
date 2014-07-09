//
// File:        KrigingModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing kriging interpolation.
//
// $Id: KrigingModel.h,v 1.2 2005/08/24 18:32:33 knap2 Exp $
//
// $Log: KrigingModel.h,v $
// Revision 1.2  2005/08/24 18:32:33  knap2
// Added data members to cache various arrays. Added getError() member
// function.
//
// Revision 1.1  2005/08/23 21:12:40  knap2
// Initial source.
//
//

#if !defined(included_krigalg_KrigingModel)
#define included_krigalg_KrigingModel

#ifndef included_iosfwd
#define included_iosfwd
#include <iosfwd>
using namespace std;
#endif

#ifndef included_krigalg_Matrix
#include <base/Matrix.h>
#endif

#ifndef included_krigalg_Value
#include <base/Value.h>
#endif

#ifndef included_krigalg_Point
#include <base/Point.h>
#endif

namespace krigalg {

  //
  // forward declarations
  //

  class KrigingModel;
  class CorrelationModel;
  class RegressionModel;

  //
  //
  //

  class KrigingModel {

    //
    // methods
    //

  public:
    //
    // construction/destruction
    //
    KrigingModel(const RegressionModel  & regressionModel,
		 const CorrelationModel & correlationModel);
    ~KrigingModel();

    //
    // meta-methods
    //

    //
    // add point/value pair to the model
    //

    void addPoint(const Point & point,
		  const Value & value);

    //
    // build the model using accumulated points
    //

    void build();

    //
    // interpolate at a point
    //

    Value interpolate(const Point & point) const;

    //
    // estimate error at a point
    //

    Value getError(const Point & point) const;

    //
    // output
    //

    friend std::ostream & operator<<(std::ostream       & outputStream,
				     const KrigingModel & krigingModel);

  private:
    //
    // copy construction/assignment
    //
    
    KrigingModel(const KrigingModel &);
    const KrigingModel & operator=(const KrigingModel &);

    bool                                  _isValid;
    std::vector<std::pair<Point, Value> > _pointData;

    Vector                                _Au;
    Vector                                _Bu;
    
    //
    // array cache
    //

    Matrix                                _matrixP;
    Matrix                                _matrixInverseR;
    Matrix                                _matrixInversePRP;

    double                                _sigmaSqr;
    
    //
    // regression/correlation models
    //

    const RegressionModel  & _regressionModel;
    const CorrelationModel & _correlationModel;

  };

}

#ifndef DEBUG_NO_INLINE
#include "KrigingModel.I"
#endif // DEBUG_NO_INLINE

#endif // included_krigalg_KrigingModel
