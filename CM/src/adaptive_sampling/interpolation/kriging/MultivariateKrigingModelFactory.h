//
// File:        MultivariateKrigingModelFactory.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Factory class for the creation of the
// Description: MultivariateKrigingModel classes.
//

#if !defined(included_krigalg_MultivariateKrigingModelFactory_h)
#define included_krigalg_MultivariateKrigingModelFactory_h
 
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

    class MultivariateKrigingModelFactory;

    //
    // local types
    //

    typedef std::shared_ptr<MultivariateKrigingModelFactory> MultivariateKrigingModelFactoryPointer;

    //
    //
    //

    /*!
     * @brief Concrete implementation of InterpolationModelFactory producing
     *        MultivariateKrigingModels
     */

    class MultivariateKrigingModelFactory : public InterpolationModelFactory {

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

      MultivariateKrigingModelFactory(const RegressionModelPointer  & regressionModel,
				      const CorrelationModelPointer & correlationModel);

      /*!
       * @brief Destructor.
       */
      
      virtual ~MultivariateKrigingModelFactory();

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

#endif // included_krigalg_MultivariateKrigingModelFactory_h
