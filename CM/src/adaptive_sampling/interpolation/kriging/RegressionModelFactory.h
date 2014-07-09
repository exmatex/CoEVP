//
// File:        RegressionModelFactory.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Factory class for the creation of the 
// Description: DerivativeRegressionModel-derived classes
//

#if !defined(included_krigalg_RegressionModelFactory_h)
#define included_krigalg_RegressionModelFactory_h

#ifndef included_krigalg_RegressionModel
#include "RegressionModel.h"
#endif

  //
  //
  //

  namespace krigalg {

    //
    // forward declarations
    //

    class RegressionModelFactory;
    
    //
    //
    //

    class RegressionModelFactory {

    public:

      //
      // types
      //

      typedef int ClassIdentifier;

      //
      // methods
      //

    public:
      //
      // construction/desctruction
      //
      RegressionModelFactory();
      ~RegressionModelFactory();
      
      //
      // meta-methods
      //
      RegressionModelPointer build(const std::string & classKey) const;
      RegressionModelPointer build(ClassIdentifier classId) const;

      static ClassIdentifier getClassId(const std::string & classKey); 
      
    private:
      //
      // copy construction and assignment operator
      //
      RegressionModelFactory(const RegressionModelFactory &);
      const RegressionModelFactory & operator=(const RegressionModelFactory &);
      
    };

}

#endif // included_krigalg_RegressionModelFactory_h
