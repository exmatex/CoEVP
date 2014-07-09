//
// File:        DerivativeRegressionModelFactory.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Factory class for the creation of the 
// Description: DerivativeRegressionModel-derived classes
//

#if !defined(included_krigalg_DerivativeRegressionModelFactory_h)
#define included_krigalg_DerivativeRegressionModelFactory_h

#ifndef included_krigalg_DerivativeRegressionModel
#include "DerivativeRegressionModel.h"
#endif

  //
  //
  //

  namespace krigalg {

    //
    // forward declarations
    //

    class DerivativeRegressionModelFactory;
    
    //
    //
    //

    class DerivativeRegressionModelFactory {

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
      DerivativeRegressionModelFactory();
      ~DerivativeRegressionModelFactory();
      
      //
      // meta-methods
      //
      DerivativeRegressionModelPointer build(const std::string & classKey) const;
      DerivativeRegressionModelPointer build(ClassIdentifier classId) const;

      static ClassIdentifier getClassId(const std::string & classKey); 
      
    private:
      //
      // copy construction and assignment operator
      //
      DerivativeRegressionModelFactory(const DerivativeRegressionModelFactory &);
      const DerivativeRegressionModelFactory & operator=(const DerivativeRegressionModelFactory &);
      
    };

}

#endif // included_krigalg_DerivativeRegressionModelFactory_h
