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
// File:        GaussianDerivativeCorrelationModel.cc
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing gaussian correlation model with derivatives.
//
// $Id: GaussianDerivativeCorrelationModel.cc,v 1.1 2005/09/01 16:14:19 knap2 Exp $
//
// $Log: GaussianDerivativeCorrelationModel.cc,v $
// Revision 1.1  2005/09/01 16:14:19  knap2
// Initial source.
//
//

#include "GaussianDerivativeCorrelationModel.h"

#include <mtl/mtl.h>

#include <cassert>
#include <cmath>
 
//
//
//

namespace krigalg {
    
  //
  // construction/destruction
  //
    
    GaussianDerivativeCorrelationModel::GaussianDerivativeCorrelationModel(const std::vector<double> & thetas)
      : DerivativeCorrelationModel(thetas)
  {

    return;

  }

  GaussianDerivativeCorrelationModel::~GaussianDerivativeCorrelationModel()
  {

    return;

  }

  //
  // get value of the correlation function between two points.
  // 
  // Z(s)={Y, D_1 Y, D_2 Y, ..., D_n Y}(s)
  //

  Matrix 
  GaussianDerivativeCorrelationModel::getValue(const Point & firstPoint,
					       const Point & secondPoint) const
  {

    //
    // compute distance between firstPoint and secondPoint
    //

    const Vector distance = firstPoint - secondPoint;


    //
    // compute L2-norm of the firstPoint-secondPoint distance
    //

    const double distanceNorm = mtl::two_norm(distance);

    //
    // allocate array 
    //

    const int arrayDimension = firstPoint.size() + 1;

    Matrix covarianceArray(arrayDimension,
			   arrayDimension);
    
    //
    // fill in values; for now the values are normalized wrt. the
    // gaussian correlation
    //

    const double theta    = _thetas.front();
    const double thetaSqr = theta*theta;
    const double covarianceScaling = exp(-theta*distanceNorm*
					 distanceNorm);

    //
    // correlation of function with itself
    //

    covarianceArray[0][0] = covarianceScaling;

    //
    // correlation of function and derivative
    //

    for (int i = 1; i < arrayDimension; ++i) {

      //
      // get covariance magnitude
      //

      const double covarianceValue = -2.0*theta*distance[i - 1];

      //
      // Cov[ D_{i-1}Y(firstPoint), Y(secondPoint) ]
      //
      
      covarianceArray[0][i] = -covarianceValue*covarianceScaling;

      //
      // Cov[ Y(firstPoint), D_{i-1} Y(secondPoint) ]
      //

      covarianceArray[i][0] = covarianceValue*covarianceScaling;

    }

    //
    // correlation of derivatives
    //

    for (int i = 1; i < arrayDimension; ++i) 
      for (int j = i; j < arrayDimension; ++j)
	if (i == j)	  
	  //
	  // self correlation
	  //
	  covarianceArray[i][i] = covarianceScaling*
	    (2.0*theta - 
	     4.0*thetaSqr*distance[i - 1]*distance[i - 1]);
	else
	  covarianceArray[i][j] = covarianceArray[j][i] =
	    covarianceScaling*(-4.0*thetaSqr*distance[i - 1]*distance[j - 1]);

    //
    // return covarianceArray
    //

    return covarianceArray;
    
  }


  //
  // get a string representation of the class name
  //

  std::string
  GaussianDerivativeCorrelationModel::getClassName() const
  {

    return std::string("krigalg::GaussianDerivativeCorrelationModel");

  }

  //
  // Database output
  //
  
  void
  GaussianDerivativeCorrelationModel::putToDatabase(toolbox::Database & db) const
  {

    //
    // store base-class data
    //
    
    DerivativeCorrelationModel::putToDatabase(db);

    //
    //
    //

    return;

  }

  //
  // Database input
  //

  void
  GaussianDerivativeCorrelationModel::getFromDatabase(toolbox::Database & db)
  {

    //
    // get base-class data
    //

    DerivativeCorrelationModel::getFromDatabase(db);

    //
    //
    //
    
    return;
    
  }

}


