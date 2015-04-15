//
// File:        MultivariateDerivativeKrigingModel.h
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

#if !defined(included_krigalg_MultivariateDerivativeKrigingModel)
#define included_krigalg_MultivariateDerivativeKrigingModel

#ifndef included_krigalg_InterpolationModel
#include <base/InterpolationModel.h>
#endif // included_krigalg_InterpolationModel

#ifndef included_krigalg_DerivativeCorrelationModel
#include "DerivativeCorrelationModel.h"
#endif // included_DerivativeCorrelationModel

#ifndef included_krigalg_DerivativeRegressionModel
#include "DerivativeRegressionModel.h"
#endif // included_DerivativeRegressionModel

#ifndef included_krigalg_Point
#include <base/Point.h>
#endif // included_krigalg_Point

#ifndef included_krigalg_Value
#include <base/Value.h>
#endif // included_krigalg_Value

#ifndef included_krigalg_Matrix
#include <base/Matrix.h>
#endif // included_krigalg_Matrix

#ifndef included_vector
#define included_vector
#include <vector>
#endif // included_vector

#ifndef included_iosfwd
#define included_iosfwd
#include <iosfwd>
#endif // included_iosfwd

//
// 
//



  //
  //
  //

  namespace krigalg {

    //
    // forward declarations
    //

    class MultivariateDerivativeKrigingModel;

    //
    // local types
    //

    typedef std::shared_ptr<MultivariateDerivativeKrigingModel> MultivariateDerivativeKrigingModelPtr;

    //
    // class definition
    //

    class MultivariateDerivativeKrigingModel : public InterpolationModel {
    
    public:
      //
      // construction
      //

      MultivariateDerivativeKrigingModel(const RegressionModelPointer  & regressionModel,
					 const CorrelationModelPointer & correlationModel);
      ~MultivariateDerivativeKrigingModel();

      //
      // features of the interpolation model.
      //

      /*!
       * Check whether the interpolation model provides error. 
       *
       * @return True if the interpolation model provides error 
       *         and false otherwise.
       */
    
      virtual bool hasError() const;

      /*!
       * Check whether the interpolation model provides gradient 
       * estimate along with the interpolated value.
       *
       * @return True if the interpolation model provides gradient
       *         and false otherwise.
       */

      virtual bool hasGradient() const;

      //
      // meta-methods
      //
    
      /*!
       * @brief Clone the object.
       *
       * @return Pointer to a copy of InterpolationModel.
       */

      virtual MultivariateDerivativeKrigingModel * clone() const;

      /*!
       * Add point and values to the model.
       * 
       * @return boolean true if the point can be added into the model and
       * false otherwise.
       *
       * @param point  Point to insert.
       * @param values Values at Point.
       */

      virtual bool addPoint(const Point              & point,
			    const std::vector<Value> & values);

      //
      // get number of points in the model
      //

      virtual int getNumberPoints() const;

      //
      // get handle to kriging model points
      //

      virtual const std::vector<Point> & getPoints() const;

      //
      // get number of values
      //

      virtual int getNumberValues() const;

      //
      // get point dimension
      //

      virtual int getPointDimension() const;

      //
      // get value dimension
      //

      virtual int getValueDimension() const;

      //
      // get regression model
      //

      RegressionModelPointer getRegressionModel() const;
    
      //
      // get correleation model
      //

      CorrelationModelPointer getCorrelationModel() const;
 
      /*!
       * Check if a model is valid.
       *
       * @return bool value; true if a model is valid and false
       * otherwise.
       */

      virtual bool isValid() const;
       
      //
      // interpolate value at a point
      //

      virtual Value interpolate(int           valueId,
				const Point & point) const;

      //
      // estimate error at a point
      //

      virtual Value getMeanSquaredError(int           valueId, 
					const Point & point) const;

      //
      // divide the current model to create two smaller models
      //

      MultivariateDerivativeKrigingModel divide();

      //
      // output
      //

      virtual void putToDatabase(toolbox::Database & db) const;
      virtual void getFromDatabase(toolbox::Database & db);
      
      virtual void pack(std::vector<double> & packedContainer) const;
      virtual void unpack(const std::vector<double> & packedContainer);

      friend std::ostream & 
	operator<<(std::ostream                             & outputStream,
		   const MultivariateDerivativeKrigingModel & krigingModel);

    protected:
      //
      // polymorphic output
      //
      virtual void print(std::ostream & outputStream) const;

    private:
      //
      // copy construction/assignment
      //
      //     MultivariateDerivativeKrigingModel(const MultivariateDerivativeKrigingModel &);
      //     const MultivariateDerivativeKrigingModel & operator=(const MultivariateDerivativeKrigingModel &);
    
      //
      // build the model using accumulated points
      //

      void build();

      //
      // data
      //

    public:

    private:
      bool                                                    _isValid;
    
      //
      // points and values
      //

      std::vector<Point>                                      _points;
      std::vector<std::vector<Value> >                        _values;
    
      //
      // value-independent data
      //

      Matrix                                                  _matrixX;
      Matrix                                                  _matrixInverseV;
      Matrix                                                  _matrixInverseXVX;
      Matrix                                                  _matrixInverseVX;

      //
      // value-dependent data
      //

      std::vector<Vector>                                    _AZ;
      std::vector<Vector>                                    _BZ;
      std::vector<double>                                    _sigmaSqr;

      //
      // regression/correlation models
      //

      RegressionModelPointer                                 _regressionModel;
      CorrelationModelPointer                                _correlationModel;

    };


}
#endif // included_krigalg_MultivariateDerivativeKrigingModel_h
