//
// File:        GaussianCorrelationModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementating gaussian correlation model.
//
// $Id: GaussianCorrelationModel.h,v 1.1 2005/08/23 21:12:40 knap2 Exp $
//
// $Log: GaussianCorrelationModel.h,v $
// Revision 1.1  2005/08/23 21:12:40  knap2
// Initial source.
//
//

#if !defined(included_krigalg_GaussianCorrelationModel)
#define included_krigalg_GaussianCorrelationModel

#ifndef included_krigalg_CorrelationModel
#include "CorrelationModel.h"
#endif

namespace krigalg {

  //
  // forward declarations
  //

  class GaussianCorrelationModel;

  //
  // class definition
  //

  class GaussianCorrelationModel : public CorrelationModel {
  
    //
    // methods
    //

  public:
    //
    // construction/destruction
    //

    GaussianCorrelationModel(const std::vector<double> & thetas);
    virtual ~GaussianCorrelationModel();

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

    GaussianCorrelationModel(const GaussianCorrelationModel &);
    const GaussianCorrelationModel & operator=(const GaussianCorrelationModel &);

  };

}

#endif // included_krigalg_GaussianCorrelationModel
