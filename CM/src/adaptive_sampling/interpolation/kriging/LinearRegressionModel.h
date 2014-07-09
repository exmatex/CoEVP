//
// File:        LinearRegressionModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing linear regression model.
//
// $Id: LinearRegressionModel.h,v 1.1 2005/08/23 21:12:40 knap2 Exp $
//
// $Log: LinearRegressionModel.h,v $
// Revision 1.1  2005/08/23 21:12:40  knap2
// Initial source.
//
//

#if !defined(included_krigalg_LinearRegressionModel)
#define included_krigalg_LinearRegressionModel

#ifndef included_krigalg_RegressionModel
#include "RegressionModel.h"
#endif

namespace krigalg {

  //
  // forward declarations
  //

  class LinearRegressionModel;

  //
  // class definition
  //

  class LinearRegressionModel : public RegressionModel {

  public:
    //
    // construction/destruction
    //

    LinearRegressionModel();
    virtual ~LinearRegressionModel();

    //
    // meta-methods
    //

    virtual Matrix    getValues(const Point         & point) const;
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
    
    LinearRegressionModel(const LinearRegressionModel &);
    const LinearRegressionModel & operator=(const LinearRegressionModel &);

  };

}

#endif // included_krigalg_LinearRegressionModel
