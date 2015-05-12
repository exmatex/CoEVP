// File:        MultivariateDerivativeKrigingModelFactory.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Factory class for the creation of the
// Description: MultivariateDerivativeKrigingModel classes.
//

#if !defined(included_krigalg_MultivariateDerivativeKrigingModelFactory_h)
#define included_krigalg_MultivariateDerivativeKrigingModelFactory_h
 
#ifndef included_krigalg_InterpolationModelFactory_h
#include <base/InterpolationModelFactory.h>
#endif 

#ifndef included_krigalg_RegressionModel
#include <kriging/RegressionModel.h>
#endif // included_krigalg_RegressionModel

#ifndef included_krigalg_CorrelationModel
#include <kriging/CorrelationModel.h>
#endif // included_krigalg_CorrelationModel

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

    class MultivariateDerivativeKrigingModelFactory;

    //
    // local types
    //

    typedef std::shared_ptr<MultivariateDerivativeKrigingModelFactory> MultivariateDerivativeKrigingModelFactoryPointer;

    //
    //
    //

    /*!
     * @brief Concrete implementation of InterpolationModelFactory producing
     *        MultivariateDerivativeKrigingModels
     */

    class MultivariateDerivativeKrigingModelFactory : public InterpolationModelFactory {

      //
      // methods
      //

    public:

      /*!
       * @brief Constructor.
       * 
       * @param regressionModel A pointer to a regression model to be used.
       * @param correlationModel A pointer to a crorrelation model to be used.
       */

      MultivariateDerivativeKrigingModelFactory(const RegressionModelPointer  & regressionModel,
						const CorrelationModelPointer & correlationModel);

      /*!
       * @brief Destructor.
       */
      
      virtual ~MultivariateDerivativeKrigingModelFactory();

      /*!
       * @brief Construct a class based on the class key.
       *
       * @return A pointer to the class instance.
       */
      
      virtual InterpolationModelPtr build() const;

      //
      // data
      //

    private:

      RegressionModelPointer  _regressionModel;
      CorrelationModelPointer _correlationModel;

    };

}

#endif // included_krigalg_MultivariateDerivativeKrigingModelFactory_h
