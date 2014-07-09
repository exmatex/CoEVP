//
// File:        DerivativeCorrelationModelFactory.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing DerivativeCorrelationModel factory
// 

#if !defined(included_krigalg_DerivativeCorrelationModelFactory_h)
#define included_krigalg_DerivativeCorrelationModelFactory_h

#ifndef included_krigalg_DerivativeCorrelationModel
#include "DerivativeCorrelationModel.h"
#endif

namespace krigalg {

    //
    // forward declarations
    //

    class DerivativeCorrelationModelFactory;

    //
    //
    //

    class DerivativeCorrelationModelFactory {

      //
      // types
      //

    public:

      typedef int ClassIdentifier;
      
      //
      // methods
      //

    public:
      //
      // construction/destruction
      //
      DerivativeCorrelationModelFactory();
      ~DerivativeCorrelationModelFactory();

      //
      // build correlation model
      //
      DerivativeCorrelationModelPointer build(const std::string & classKey) const;
      DerivativeCorrelationModelPointer build(const std::string & classKey,
					      const std::vector<double> & thetas) const;
      
      DerivativeCorrelationModelPointer build(ClassIdentifier classId) const;
      DerivativeCorrelationModelPointer build(ClassIdentifier classId,
					      const std::vector<double> & thetas) const;

      static ClassIdentifier getClassId(const std::string & classKey);

    private:
      //
      // copy construction and assignment operator
      //
      DerivativeCorrelationModelFactory(const DerivativeCorrelationModelFactory &);
      const DerivativeCorrelationModelFactory & operator=(const DerivativeCorrelationModelFactory &);

    };

}

#endif // included_krigalg_DerivativeCorrelationModelFactory_h
