//
// File:        DerivativeKrigingModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing kriging model including derivative 
// Description: information.
//
// $Id$
//
// $Log$
//

#if !defined(included_krigalg_DerivativeKrigingModel)
#define included_krigalg_DerivativeKrigingModel

#ifndef included_iosfwd
#define included_iosfwd
#include <iosfwd>
using namespace std;
#endif

#ifndef included_krigalg_DerivativeCorrelationModel
#include "DerivativeCorrelationModel.h"
#endif // included_DerivativeCorrelationModel

#ifndef included_krigalg_DerivativeRegressionModel
#include "DerivativeRegressionModel.h"
#endif // included_krigalg_DerivativeRegressionModel

#ifndef included_krigalg_Point
#include <base/Point.h>
#endif

#ifndef included_krigalg_Value
#include <base/Value.h>
#endif

#ifndef included_krigalg_Matrix
#include <base/Matrix.h>
#endif

namespace krigalg {

  //
  // class definition
  //

  class DerivativeKrigingModel {
    
    //
    // methods
    //
    
  public:
    //
    // construction/destruction
    //
    
    DerivativeKrigingModel(const DerivativeRegressionModelPointer  & regressionModel,
			   const DerivativeCorrelationModelPointer & correlationModel);
    ~DerivativeKrigingModel();

    //
    // meta-methods
    //

    //
    // add point/value pair to the model
    //

    void addPoint(const Point & point,
		  const Value & value);

    //
    // get number of points in the model
    //

    int getNumberPoints() const;

    //
    // get kriging model points
    //

    std::vector<Point> getPoints() const;

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

    Value getMeanSquaredError(const Point & point) const;

    //
    // divide the current model to create two smaller models
    //

    DerivativeKrigingModel divide();

    //
    // output
    //

    friend std::ostream & operator<<(std::ostream                 & outputStream,
				     const DerivativeKrigingModel & krigingModel);

  private:
    //
    // copy construction/assignment
    //
    
    // DerivativeKrigingModel(const DerivativeKrigingModel &);
    //const DerivativeKrigingModel & operator=(const DerivativeKrigingModel &);

  private:
    bool                                  _isValid;
    std::vector<std::pair<Point, Value> > _pointData;
    
    Vector                                _AZ;
    Vector                                _BZ;
    
    //
    // array cache
    //

    Matrix                                _matrixX;
    Matrix                                _matrixInverseV;
    Matrix                                _matrixInverseXVX;

    double                                _sigmaSqr;
    
    //
    // regression/correlation models
    //

    DerivativeRegressionModelPointer  _regressionModel;
    DerivativeCorrelationModelPointer _correlationModel;

  };

}

#endif // included_krigalg_DerivativeKrigingModel
