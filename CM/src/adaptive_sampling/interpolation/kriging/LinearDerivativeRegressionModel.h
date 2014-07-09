//
// File:        LinearDerivativeRegressionModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing linear regression model with derivatives.
//
// $Id$
//
// $Log$
//

#if !defined(included_krigalg_LinearDerivativeRegressionModel)
#define included_krigalg_LinearDerivativeRegressionModel

#ifndef included_krigalg_DerivativeRegressionModel
#include "DerivativeRegressionModel.h"
#endif

namespace krigalg {

  //
  // forward declarations
  //

  class LinearDerivativeRegressionModel;

  //
  // class definition
  //

  class LinearDerivativeRegressionModel : public DerivativeRegressionModel {

    //
    // methods
    //

  public:
    //
    // construction/destruction
    //
    LinearDerivativeRegressionModel();
    virtual ~LinearDerivativeRegressionModel();
    
    //
    // meta-methods
    //
    
    virtual Matrix    getValues(const Point & point) const;
    virtual Dimension getDimension(const Point & point) const;

    //
    // Database input/output
    //

    virtual std::string getClassName() const;
    virtual void putToDatabase(toolbox::Database & db) const;
    virtual void getFromDatabase(toolbox::Database & db);

  private:
    //
    // copy construction/assignment
    //
    LinearDerivativeRegressionModel(const LinearDerivativeRegressionModel &);
    const LinearDerivativeRegressionModel & operator=(const LinearDerivativeRegressionModel &);

  };

}

#endif // included_krigalg_LinearDerivativeRegressionModel
