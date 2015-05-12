//
// File:        CorrelationModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Abstract base class for stochastic correlation model
//
// $Id: CorrelationModel.h,v 1.1 2005/08/23 21:12:40 knap2 Exp $
//
// $Log: CorrelationModel.h,v $
// Revision 1.1  2005/08/23 21:12:40  knap2
// Initial source.
//
//

#if !defined(included_krigalg_CorrelationModel)
#define included_krigalg_CorrelationModel

#ifndef included_vector
#define included_vector
#include <vector>
using namespace std;
#endif

#include <memory>

#ifndef included_krigalg_Point
#include <base/Point.h>
#endif

#ifndef included_krigalg_Matrix
#include <base/Matrix.h>
#endif

  //
  // forward declarations
  //

  namespace toolbox {

    class Database;
    
  }

  //
  //
  //

  namespace krigalg {

  //
  // forward declarations
  //

  class CorrelationModel;
  typedef std::shared_ptr<CorrelationModel> CorrelationModelPointer;

  //
  // class definition
  //

  class CorrelationModel {

    //
    // methods
    //

  public:
    //
    // construction/destruction
    //

    CorrelationModel(const std::vector<double> & thetas);
    virtual ~CorrelationModel() = 0;

    //
    // meta-methods
    //

    virtual Matrix getValue(const Point & firstPoint,
			    const Point & secondPoint) const = 0;
    void getThetas(std::vector<double> & thetas) const;
    void setThetas(const std::vector<double> & thetas);

    //
    // Database input/output
    //

    virtual std::string getClassName() const = 0;
    virtual void putToDatabase(toolbox::Database & db) const;
    virtual void getFromDatabase(toolbox::Database & db);

  protected:
    std::vector<double> _thetas;

  private:
    //
    // copy construction/assignment
    //
    
    CorrelationModel(const CorrelationModel &);
    const CorrelationModel & operator=(const CorrelationModel &);

  };

}

#endif // krigalg_CorrelationModel
