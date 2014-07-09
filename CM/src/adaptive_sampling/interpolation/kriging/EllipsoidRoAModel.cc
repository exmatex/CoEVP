// DO-NOT-DELETE revisionify.begin() 
/*

                            Copyright (c) 2014.
               Lawrence Livermore National Security, LLC.
         Produced at the Lawrence Livermore National Laboratory
                             LLNL-CODE-656392.
                           All rights reserved.

This file is part of CoEVP, Version 1.0. Please also read this link -- http://www.opensource.org/licenses/index.php

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.

*/
// DO-NOT-DELETE revisionify.end() 
//
// File:        EllipsoidRoAModel.cc
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

#include "EllipsoidRoAModel.h"

// #include "SecondMoment.h"

#include <mtl/mtl.h>
#include <mtl/mtl2lapack.h>

#include <algorithm>
#include <cassert>
#include <iostream>
// #include <iterator>
// #include <limits>

//
//
// 

namespace ellalg {

     //
     //
     //

     bool EllipsoidRoAModel::_haveComputedPrivateStuff = false;

     // for ellipsoid growth control:
     double EllipsoidRoAModel::_grH = 1.2;
     double EllipsoidRoAModel::_grM = 64.0;
     double EllipsoidRoAModel::_grK = -1.; // computed from _grH and _grM
     double EllipsoidRoAModel::_grMInv = -1.; // computed from _grH and _grM
     bool   EllipsoidRoAModel::_grExtend = false; 

     // for ellipsoid shape control:
     double EllipsoidRoAModel::_shM = -1.;
     bool   EllipsoidRoAModel::_shZero = false;
     bool   EllipsoidRoAModel::_shInf = false;

     // for extra interpolation control
     double EllipsoidRoAModel::_interpPDistMax = 1.0;

    // for initial A matrix size
    double EllipsoidRoAModel::_epsA = 1.0e-5;
     

    //
    // construction/destruction
    //
    EllipsoidRoAModel::EllipsoidRoAModel(
					 const krigalg::Matrix & GJ,
					 const krigalg::Point & point, 
					 const krigalg::Value & value
					 )
      : _GJ(GJ),
	_point(point), 
	_value(value),
	_A(value.size(), value.size()) // initializes entries to zero
    {

      if (!(EllipsoidRoAModel::_haveComputedPrivateStuff)) {
	computePrivateStuff();
      }

      mtl::set_diagonal(_A, 1./EllipsoidRoAModel::_epsA);

      return;

    }

    EllipsoidRoAModel::~EllipsoidRoAModel()
    {

      return;

    }

    //
    // set precomputed private stuff
    //
    void
    EllipsoidRoAModel::computePrivateStuff(){

       if (EllipsoidRoAModel::_grH < 0.) {
	  EllipsoidRoAModel::_grExtend = false;
	  _haveComputedPrivateStuff = true;
	  return;
       }

       assert(EllipsoidRoAModel::_grH >= 1.);
       
       EllipsoidRoAModel::_grK = pow(EllipsoidRoAModel::_grH, EllipsoidRoAModel::_grM) -1.; 
       EllipsoidRoAModel::_grMInv = 1. / EllipsoidRoAModel::_grM;

       EllipsoidRoAModel::_grExtend = true;

       if (EllipsoidRoAModel::_shM < 0.) {
	  EllipsoidRoAModel::_shInf  = true;
	  EllipsoidRoAModel::_shZero = false;
       } else if (EllipsoidRoAModel::_shM < 1e-8) {
	  EllipsoidRoAModel::_shInf = false;
	  EllipsoidRoAModel::_shZero = true;
       } 

       _haveComputedPrivateStuff = true;

      return;
    }

    //
    // set parameters
    //
    void
    EllipsoidRoAModel::setParams(double shM, double grH, double grM, double interpPDistMax, double epsA) 
    {
       EllipsoidRoAModel::_shM = shM;
       EllipsoidRoAModel::_grH = grH;
       EllipsoidRoAModel::_grM = grM;
       EllipsoidRoAModel::_interpPDistMax = interpPDistMax;
       EllipsoidRoAModel::_epsA = epsA;
       computePrivateStuff();
    }

    //
    // get point dimension
    //
  
    int 
    EllipsoidRoAModel::getPointDimension() const
    {

      return _point.size();

    }
  
    //
    // get value dimension
    //
  
    int 
    EllipsoidRoAModel::getValueDimension() const
    {

      return _value.size();

    }

    krigalg::Value
    EllipsoidRoAModel::interpolate(const krigalg::Point & point) const
    {

      const krigalg::Vector inputDiff = point - _point;

      const krigalg::Vector estValueDiff = mult(_GJ, inputDiff);
      
      const krigalg::Vector interpolatedValueVector = _value + estValueDiff;
      return interpolatedValueVector;

    }

    krigalg::Value
    EllipsoidRoAModel::interpolate(const krigalg::Value & valueDiff) const
    {

      const krigalg::Vector interpolatedValueVector = _value + valueDiff;
      return interpolatedValueVector;

    }

    double
    EllipsoidRoAModel::determineGrowthFactor(const double errorRatio) const {
      assert(errorRatio <= 1.);
      assert(errorRatio >= 0.);

      double growthFactor;
	  
      if (_grExtend) {
	growthFactor = _grH * pow( (1./(1.+_grK * errorRatio)), _grMInv); 
      } else {
	growthFactor = 1.;
      }
	  
      return growthFactor;
    }

    void
    EllipsoidRoAModel::doEllipsoidGrowth(
		      const krigalg::Value & valueDiff, // to be enclosed
		      const krigalg::Point & pointNew, // point at which new evaluation was done
		      const double errorRatio, // ratio of actual to acceptable error
		      bool & shifted, // whether or not shifted center
		      double & shiftFactor // amount by which shifted
		      ){

	  const int valueDimension = getValueDimension();

	  krigalg::Vector nMap = mult(_A, valueDiff);

	  const double dSq = mtl::dot( nMap, nMap );

	  if (dSq < 1.) {
	     // may be possible under certain circumstances -- checking for growth of a whole bunch
	     // of models, using a point that may be inside some of them already
	     shifted = false ;
	     shiftFactor = 0. ;
	     return;
	  }
	  
	  double d = sqrt(dSq); // = mtl::two_norm(nMap);
	  mtl::scale(nMap, 1./d);

	  const double growthFactor = determineGrowthFactor(errorRatio);
	  // with GrowthFactor > 1, can grow the ellipsoid a little more than is strictly needed to touch the new value
	  d *= growthFactor ;

	  if (_shInf) {
	     // in limit m->inf, no lateral growth

	     const double factor = (1. / d - 1.);

	     for (int jA = 0; jA < valueDimension; jA++){
	       const double tempknA = factor * mtl::dot( nMap, columns(_A)[jA] ) ;
	       // factor * SUM(n_map(:)*_A(:,j_A)); // from f90 notation

	       // columns(_A)[jA] += nMap * tempknA;
	       mtl::add(columns(_A)[jA],
			mtl::scaled(nMap, tempknA),
			columns(_A)[jA]);
	       // columns(_A)[jA] += nMap * tempknA; // does not work
	       // A_matx(:,j_A) = A_matx(:,j_A) + n_map(:)*temp_knA // from f90 notation
	     }

	     shifted = false ;
	     shiftFactor = 0. ;

	  } else {

	     double a, b;

	     if (_shZero) {
		// in limit m->0, no backside growth
		a = 0.5 * (d + 1.);
		b = sqrt(a);
	     } else {
		const double dd = d*d;
		a = ((2.+_shM)*d + sqrt(4.*(1.+_shM) + dd*_shM*_shM ))/(4. + 2.*_shM);
		const double twoad = 2.*a*d;
		const double factor_b = 1. - dd + twoad;
		const double temp = factor_b*factor_b - 4.*a*a;
		if (temp <= 0.) {
		   b = sqrt( factor_b * 0.5) ;
		} else {
		   b = sqrt(( factor_b - sqrt(temp) )*0.5);
		}
	     }

	     const double capX = d - a;

	     const double bInv =   1. / b;
	     const double factor   = ((1. / a) - bInv);
	     for (int jA = 0; jA < valueDimension; jA++){
	       const double tempknA = factor * mtl::dot( nMap, columns(_A)[jA] ) ;
	       // factor * SUM(n_map(:)*_A(:,j_A)); // from f90 notation
	       // columns(_A)[jA] *= bInv;
	       mtl::scale(columns(_A)[jA], bInv) ;
	       mtl::add(columns(_A)[jA],
			mtl::scaled(nMap, tempknA),
			columns(_A)[jA]);
	       // columns(_A)[jA] += nMap * tempknA;  // does not work
	       // _A(:,j_A) = bInv * _A(:,j_A) + n_map(:)*tempknA; // from f90 notation
	     }

	     shifted    = true ;
	     shiftFactor = capX/d;

	  }

	  if (shifted) {

	    // shift intput and output (point and value) vectors;
	    // can do this because using a linear approximation;
	    // all the mappings from input vector (point) through to get shiftFactor
	    // are linear, so shiftFactor can be applied to estValueDiff and to
	    // inputDiff;

	    // shift point toward point at which new evaluation was done
	    _point *= (1. - shiftFactor);
	    _point += pointNew * shiftFactor;

	    // shift reference value by the amount dictated by GJ
	    _value += valueDiff * shiftFactor;

	  }

	  return ;
       }

       void
       EllipsoidRoAModel::testInterpG(
		   const krigalg::Point & inputPoint,
		   bool & canInterp,
		   bool & hitLimitIDist,
		   krigalg::Value & estOutputDiff
		   ){

	  canInterp = false ;
	  hitLimitIDist = false ;

	  const krigalg::Vector inputDiff      =  inputPoint - _point;
	  estOutputDiff  = mult(_GJ, inputDiff); // _GJ is in normalized spaces
	  const double inputPDist = mtl::two_norm(inputDiff); // already normalized

	  if (inputPDist <= _interpPDistMax) {
	    // check ellipsoid;
	    // estOutputDiff already normalized;
	    const krigalg::Vector outputInterpDistA = mult(_A, estOutputDiff);
	     const double oDistA = mtl::two_norm(outputInterpDistA);
	     if (oDistA <= 1.) {
		canInterp = true;
	     }
	  } else {
	     hitLimitIDist = true;
	     canInterp = false;
	  }

	  return ;
       }


    //
    // get current ellipsoid center
    //

    krigalg::Point 
    EllipsoidRoAModel::getCenter() const
    {

      return _point;

    }

    //
    // get GJ
    //

     void
     EllipsoidRoAModel::getGJ(double * gradient) const
     {

	const krigalg::Matrix gjTranspose = krigalg::transpose(_GJ);

	const double * gjData = gjTranspose.data();
	std::copy(&(gjData[0]), &(gjData[_value.size()*_point.size()]), gradient);

	return;

    }

    //
    // output object data to database; store only point-value pairs and
    // rely in the getFromDatabase() to properly build the model
    //

    void
    EllipsoidRoAModel::putToDatabase(toolbox::Database& db) const
    {

      assert(false); // Needs to be implemented
      return;

    }

    //
    // fill in the object from Database
    //

    void
    EllipsoidRoAModel::getFromDatabase(toolbox::Database& db)
    {

      assert(false); // Needs to be implemented
      return;

    }
      

    //
    // output operator
    //

    std::ostream &
    operator<<(std::ostream                             & outputStream,
	       const EllipsoidRoAModel & ellipsoidModel)
    {

      outputStream << "point: " << std::endl;
      outputStream << ellipsoidModel._point << std::endl;

      outputStream << "value: " << std::endl;
      outputStream << ellipsoidModel._value << std::endl;

      outputStream << "GJ: " << std::endl;
      //outputStream << ellipsoidModel._GJ;
      krigalg::operator<<(outputStream, ellipsoidModel._GJ);
      outputStream << std::endl;

      outputStream << "A: " << std::endl;
      // outputStream << ellipsoidModel._A;
      krigalg::operator<<(outputStream, ellipsoidModel._A);
      outputStream << std::endl;

      return outputStream;

    }

}


