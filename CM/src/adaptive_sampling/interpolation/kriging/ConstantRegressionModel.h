//
// File:        ConstantRegressionModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Derived class for the constant regression model used in 
// Description: Kriging computations.
//
// $Id: ConstantRegressionModel.h,v 1.1 2005/08/23 21:12:40 knap2 Exp $
// 
// $Log: ConstantRegressionModel.h,v $
// Revision 1.1  2005/08/23 21:12:40  knap2
// Initial source.
//
//

#if !defined(included_krigalg_ConstantRegressionModel)
#define included_krigalg_ConstantRegressionModel

#ifndef included_krigalg_RegressionModel
#include "RegressionModel.h"
#endif

namespace krigalg {

  //
  // forward declarations
  //

  class ConstantRegressionModel;

  //
  // class definition
  //

  class ConstantRegressionModel : public RegressionModel {

  public:
    //
    // construction/destruction
    //
    
    ConstantRegressionModel();
    virtual ~ConstantRegressionModel();

    //
    // meta-methods
    //

    virtual Matrix    getValues(const Point & point) const;
    virtual Dimension getDimension(const Point & point) const;

    //
    // Database input/output
    //

    virtual std::string getClassName() const = 0;
    virtual void putToDatabase(toolbox::Database & db) const;
    virtual void getFromDatabase(toolbox::Database & db);

  private:
    //
    // copy construction/assignment
    //
 
    ConstantRegressionModel(const ConstantRegressionModel &);
    const ConstantRegressionModel & operator=(const ConstantRegressionModel &);

  };

}

#endif // included_krigalg_ConstantRegressionModel
