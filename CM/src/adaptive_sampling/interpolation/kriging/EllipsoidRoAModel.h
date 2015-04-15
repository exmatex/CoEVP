//
// File:        EllipsoidRoAModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing ellipsoid RoA model
//
// $Id$
//
// $Log$
//

#if !defined(included_ellalg_EllipsoidRoAModel)
#define included_ellalg_EllipsoidRoAModel

#ifndef included_krigalg_Point
#include <base/Point.h>
#endif // included_krigalg_Point

#ifndef included_krigalg_Value
#include <base/Value.h>
#endif // included_krigalg_Value

#ifndef included_krigalg_Matrix
#include <base/Matrix.h>
#endif // included_krigalg_Matrix

#include <memory>

#ifndef included_vector
#define included_vector
#include <vector>
#endif // included_vector

#ifndef included_iosfwd
#define included_iosfwd
#include <iosfwd>
#endif // included_iosfwd

//
// 
//

  //
  // forward declarations
  //

  namespace toolbox {

    class Database;
    
  }

  namespace ellalg {

    //
    // forward declarations
    //

    class EllipsoidRoAModel;

    //
    // local types
    //

    typedef std::shared_ptr<EllipsoidRoAModel> EllipsoidRoAModelPtr;

    //
    // class definition
    //

    class EllipsoidRoAModel {
    
    public:
      //
      // construction
      //

      EllipsoidRoAModel(
			const krigalg::Matrix & GJ,
			const krigalg::Point & point, 
			const krigalg::Value & value
			);
      ~EllipsoidRoAModel();

      //
      // meta-methods
      //
    
      //
      // set parameters
      //
      static void setParams(double shM, double grH, double grM, double interpPDistMax, double epsA) ;

      //
      // get point dimension
      //

      int getPointDimension() const;

      //
      // get value dimension
      //

      int getValueDimension() const;

      //
      // interpolate value at a point
      //

      krigalg::Value interpolate(const krigalg::Point & point) const;
      krigalg::Value interpolate(const krigalg::Value & valueDiff ) const;

      //
      // do ellipsoid growth
      //
      void doEllipsoidGrowth(
			     const krigalg::Value & valueDiff,
			     const krigalg::Point & pointNew,
			     const double errorRatio,
			     bool & shifted,
			     double & shiftFactor);

      //
      // test to see if can interpolate
      //
      void testInterpG(
		       const krigalg::Point & inputPoint,
		       bool & canInterp,
		       bool & hitLimitIDist,
		       krigalg::Value & estOutputDiff);

      //
      // get current ellipsoid center
      //

      krigalg::Point getCenter() const;

      //
      // get GJ
      //

      void getGJ(double * gradient) const;

      //
      // output
      //

      void putToDatabase(toolbox::Database & db) const;
      void getFromDatabase(toolbox::Database & db);

      friend std::ostream & 
	operator<<(std::ostream                             & outputStream,
		   const EllipsoidRoAModel & ellipsoidModel);

    private:
      //
      // copy construction/assignment
      //
      //     EllipsoidRoAModel(const EllipsoidRoAModel &);
      //     const EllipsoidRoAModel & operator=(const EllipsoidRoAModel &);
    

      //
      // set precomputed private stuff
      //
      static void computePrivateStuff();

      double determineGrowthFactor(const double errorRatio) const;

      

    public:

    private:
    
      krigalg::Point        _point; // reference input
      krigalg::Value        _value; // reference output

      
      krigalg::Matrix   _GJ; // generalized Jacobian
      krigalg::Matrix   _A; // used in output distance measure

      static bool _haveComputedPrivateStuff;

      // for ellipsoid growth control:
      static double _grH, _grK, _grM, _grMInv;
      static bool   _grExtend; 

      // for ellipsoid shape control:
      static double _shM; // NOT dimension of space; this is exponent on minor axis in volume expression
      static bool   _shZero, _shInf;

      // for extra interpolation control
      static double _interpPDistMax;

      // for initial A matrix size
      static double _epsA;


    };

}
#endif // included_ellalg_EllipsoidRoAModel_h
