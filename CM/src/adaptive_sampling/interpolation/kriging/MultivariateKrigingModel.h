//
// File:        MultivariateKrigingModel.h
// Package:     kriging algorithm
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: Class implementing multivariate kriging model.
//

#if !defined(included_krigalg_MultivariateKrigingModel)
#define included_krigalg_MultivariateKrigingModel

#ifndef included_krigalg_InterpolationModel
#include <base/InterpolationModel.h>
#endif // included_krigalg_InterpolationModel

#ifndef included_krigalg_Matrix
#include <base/Matrix.h>
#endif // included_krigalg_Matrix

#ifndef included_krigalg_CorrelationModel
#include "CorrelationModel.h"
#endif // included_CorrelationModel

#ifndef included_krigalg_RegressionModel
#include "RegressionModel.h"
#endif // included_RegressionModel

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

    class MultivariateKrigingModel;

    //
    // local types
    //

    typedef std::shared_ptr<MultivariateKrigingModel> MultivariateKrigingModelPtr;

    //
    // class definition
    //

    class MultivariateKrigingModel : public InterpolationModel {

    public:
      //
      // construction
      //

      MultivariateKrigingModel(const RegressionModelPointer  & regressionModel,
			       const CorrelationModelPointer & correlationModel);
      virtual ~MultivariateKrigingModel();

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
      
      /*!
       * @brief Clone the object.
       *
       * @return Pointer to a copy of InterpolationModel.
       */
      
      virtual MultivariateKrigingModel * clone() const;

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
      
      /*!
       * Get number of points in the model.
       *
       * @return integer number of points currently in the model.
       */
      
      virtual int getNumberPoints() const;

      /*!
       * Get handle to model points.
       *
       * @return const reference to an STL vector holding the point data.
       */

      virtual const std::vector<Point> & getPoints() const;

      /*!
       * Get number of values in the model.
       *
       * @return integer number of values.
       */

      virtual int getNumberValues() const;

      /*!
       * Get point dimension.
       *
       * @return integer point dimension.
       */

      virtual int getPointDimension() const;

      /*!
       * Get value dimension.
       *
       * @return integer value dimension.
       */

      virtual int getValueDimension() const;

      /*!
       * Get regression model.
       *
       * @return object RegressionModelPointer used by the model. 
       */

      RegressionModelPointer getRegressionModel() const;
    
      /*!
       * Get correleation model.
       *
       * @return object CorrelationModelPointer used by the model.
       */

      CorrelationModelPointer getCorrelationModel() const;
        
      /*!
       * Check if a model is valid.
       *
       * @return bool value; true if a model is valid and false
       * otherwise.
       */

      virtual bool isValid() const;

      /*!
       * Interpolate value at a point.
       *
       * @param valueId index of the value to interpolate.
       * @param point reference to a point to interpolate at.
       * 
       * @return value at the point.
       */

      virtual Value interpolate(int           valueId,
				const Point & point) const;

      /*!
       * Estimate error at a point.
       *
       * @param valueId index of the value to interpolate.
       * @param point reference to a point to interpolate at.
       * 
       * @return mean squared error at the point.
       */

      virtual Value getMeanSquaredError(int           valueId, 
					const Point & point) const;

      /*!
       * Subdivide the current model to create two smaller models.
       *
       * @return MultivariateKrigingModel object.
       */

      MultivariateKrigingModel divide();

      //
      // output
      //

      virtual void putToDatabase(toolbox::Database & db) const;
      virtual void getFromDatabase(toolbox::Database & db);
      
      virtual void pack(std::vector<double> & packedContainer) const;
      virtual void unpack(const std::vector<double> & packedContainer);

      friend std::ostream & 
	operator<<(std::ostream                   & outputStream,
		   const MultivariateKrigingModel & krigingModel);

    protected:
      //
      // polymorphic output
      //
      virtual void print(std::ostream & outputStream) const;

    private:
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

#endif // included_krigalg_MultivariateKrigingModel
