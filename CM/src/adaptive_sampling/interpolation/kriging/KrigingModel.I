/* DO-NOT-DELETE revisionify.begin() */
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
/* DO-NOT-DELETE revisionify.end() */
//
// $Id: KrigingModel.i.cc,v 1.2 2005/08/24 18:33:37 knap2 Exp $
//
// $Log: KrigingModel.i.cc,v $
// Revision 1.2  2005/08/24 18:33:37  knap2
// Added computeTranspose() function. Initial implementation of
// KrigingModel::getError().
//
// Revision 1.1  2005/08/23 21:12:40  knap2
// Initial source.
//
//

#ifndef included_krigalg_CorrelationModel
#include "CorrelationModel.h"
#endif

#ifndef included_krigalg_RegressionModel
#include "RegressionModel.h"
#endif

#ifndef included_mtl_mtl
#define included_mtl_mtl
#include <mtl/mtl.h>
#endif

#ifndef included_cassert
#define included_cassert
#include <cassert>
#endif

#ifndef included_numeric
#define included_numeric
#include <numeric>
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif // DEBUG_NO_INLINE

namespace krigalg {
  
    //
    // interpolate
    //

    inline Value 
    KrigingModel::interpolate(const Point & point) const
    {

      //
      // firewalls
      //

      assert(_isValid == true);

      //
      //
      //

      Value interpolatedValue(1);

      //
      // obtain regression model values at point
      //

      const Matrix regressionModelValuesTemp = 
	_regressionModel.getValues(point);

      assert(regressionModelValuesTemp.ncols() == 1);

      Vector regressionModelValues(regressionModelValuesTemp.nrows());

      for (int i = 0; i < regressionModelValuesTemp.nrows(); ++i)
	regressionModelValues[i] = regressionModelValuesTemp[i][0];

      //
      // obtain correlation model values between point and points
      // stored in _pointData
      //

      std::vector<double> correlationModelValues;
      correlationModelValues.reserve(_Bu.size());

      std::vector<std::pair<Point, Value> >::const_iterator pointDataIter;
      std::vector<std::pair<Point, Value> >::const_iterator pointDataEnd =
	_pointData.end();

      for (pointDataIter  = _pointData.begin();
	   pointDataIter != pointDataEnd;
	   ++pointDataIter)
	correlationModelValues.push_back((_correlationModel.getValue((*pointDataIter).first,
								     point))[0][0]);

      //
      //
      //

      assert(regressionModelValues.size()  == _Au.size());
      assert(correlationModelValues.size() == _Bu.size());

      interpolatedValue[0] = 
	std::inner_product(regressionModelValues.begin(),
			   regressionModelValues.end(),
			   _Au.begin(),
			   0.0) +
	std::inner_product(correlationModelValues.begin(),
			   correlationModelValues.end(),
			   _Bu.begin(),
			   0.0);

      //
      //
      //

      return interpolatedValue;

    }
  
    //
    // get (estimated) interpolation error
    //

    inline Value
    KrigingModel::getError(const Point & point) const
    {

      //
      // firewalls
      //

      assert(_isValid == true);
    
      //
      // obtain correlation model values between point and points
      // stored in _pointData
      // 

      Vector r(_pointData.size());

      std::vector<std::pair<Point, Value> >::const_iterator pointDataIter;
      std::vector<std::pair<Point, Value> >::const_iterator pointDataEnd =
	_pointData.end();
      int valuesOffset = 0;

      for (pointDataIter  = _pointData.begin();
	   pointDataIter != pointDataEnd;
	   ++pointDataIter, ++valuesOffset)
	r[valuesOffset] = (_correlationModel.getValue((*pointDataIter).first,
						      point))[0][0];

      //
      // obtain regression model values at point
      //

      const Matrix regressionModelValuesArrayTemp = 
	_regressionModel.getValues(point);

      assert(regressionModelValuesArrayTemp.ncols() == 1);

      Vector regressionModelValuesArray(regressionModelValuesArrayTemp.nrows()); 
      
      for (int i = 0; i < regressionModelValuesArrayTemp.nrows(); ++i)
	regressionModelValuesArray[i] = regressionModelValuesArrayTemp[i][0];
      

      //
      // compute u
      // 

      const Matrix transposeMatrixP = transpose(_matrixP);
      const Matrix tmpMatrix = mult(transposeMatrixP, 
				    _matrixInverseR);
      Vector tmpVector = mult(tmpMatrix,
			      r);

      const Vector u(tmpVector.size());
      mtl::copy(tmpVector - regressionModelValuesArray, u);

      //
      // compute u.Inverse[Transpose[P].Inverse[R].P].u 
      // 

      tmpVector = mult(_matrixInversePRP, u);
      const double part1 = mtl::dot(u, tmpVector);

      //
      // compute r.Inverse[R].r
      //
    
      tmpVector = mult(_matrixInverseR, r);
      const double part2 = mtl::dot(r, tmpVector);

      const double errorPart = part1 - part2;

      //
      // instantiate return value
      //
    
      Value error(1);

      //
      // compute error estimate
      //

      error[0] = _sigmaSqr*(1.0 + errorPart);
     
      //
      //
      //
    
      return error;

    }

}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif // DEBUG_NO_INLINE


