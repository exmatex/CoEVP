//
// File:        GaussianDerivativeCorrelationModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing gaussian correlation model with derivatives.
//
// $Id: GaussianDerivativeCorrelationModel.h,v 1.1 2005/09/01 16:14:46 knap2 Exp $
//
// $Log: GaussianDerivativeCorrelationModel.h,v $
// Revision 1.1  2005/09/01 16:14:46  knap2
// Initial source.
//
//

#if !defined(included_krigalg_GaussianDerivativeCorrelationModel)
#define included_krigalg_GaussianDerivativeCorrelationModel

#ifndef included_krigalg_DerivativeCorrelationModel
#include "DerivativeCorrelationModel.h"
#endif

namespace krigalg {

  //
  // forward declarations
  //

  class GaussianDerivativeCorrelationModel;

  //
  // class definition
  //

  class GaussianDerivativeCorrelationModel : public DerivativeCorrelationModel 
  {
  
    //
    // methods
    //

  public:
    //
    // construction/destruction
    //

    GaussianDerivativeCorrelationModel(const std::vector<double> & thetas);
    virtual ~GaussianDerivativeCorrelationModel();

    //
    // meta-methods
    //

    virtual Matrix getValue(const Point & firstPoint,
			    const Point & secondPoint) const;

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

    GaussianDerivativeCorrelationModel(const GaussianDerivativeCorrelationModel &);
    const GaussianDerivativeCorrelationModel & operator=(const GaussianDerivativeCorrelationModel &);

  };

}

#endif // included_krigalg_GaussianDerivativeCorrelationModel
