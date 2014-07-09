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
// File:        DerivativeKrigingModel.cc
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing kriging model including derivative 
// Description: information.
//
// $Id$
//
// $Log$
//

#include "DerivativeKrigingModel.h"

#include "DerivativeCorrelationModel.h"
#include "DerivativeRegressionModel.h"
#include "SecondMoment.h"

#include <mtl/mtl.h>
#include <mtl/mtl2lapack.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <limits>
#include <iostream>

//
//
//

namespace krigalg {

    //
    //
    // 

    namespace {
    
      //
      // reload values at points into a matrix
      //

      Vector 
      reloadValues(const std::vector<std::pair<Point, Value> > & _pointData)
      {
      
	//
	// firewalls
	//

	assert(_pointData.empty() == false);

	//
	// get point/value dimension
	//

	const int pointDimension = _pointData.front().first.size();
	const int valueDimension = _pointData.front().second.size();     
	const int numberPoints = _pointData.size();

	//
	// firewalls
	//

	assert(pointDimension + 1 == valueDimension);

	//
	// instantiate
	//

	Vector Z(valueDimension*numberPoints);

	//
	// insert values data into Z
	//

	std::vector<std::pair<Point, Value> >::const_iterator pointDataIter;
	std::vector<std::pair<Point, Value> >::const_iterator pointDataEnd = 
	  _pointData.end();
	int currentPoint = 0;

	for (pointDataIter  = _pointData.begin();
	     pointDataIter != pointDataEnd;
	     ++pointDataIter, ++currentPoint) {
	
	  //
	  // get handle to Value
	  //

	  const Value & value = (*pointDataIter).second;
	  assert(value.size() == valueDimension);

	  //
	  // copy Value
	  //

	  for (Value::size_type i = 0; i < value.size(); ++i)
	    Z[currentPoint + i*numberPoints] = value[i];

	}

	return Z;

      }

      //
      // construct correlation matrix for all points in the model
      //   

      Matrix 
      createCorrelationMatrix(const std::vector<std::pair<Point, Value> > & _pointData,
			      DerivativeCorrelationModelPointer & _correlationModel)
      {

	//
	// firewalls
	//

	assert(_pointData.empty() == false);

	//
	// get dimensions
	//
      
	const int pointDimension = _pointData.front().first.size();
	const int valueDimension = _pointData.front().second.size();
	const int numberPoints   = _pointData.size();

	//
	// instantiate matrix
	//

	Matrix V(valueDimension*numberPoints,
		 valueDimension*numberPoints);

	//
	// iterate over all points
	//

	for (int i = 0; i < numberPoints; ++i) {

	  //
	  // get i-point handle 
	  //

	  const Point & iPoint = _pointData[i].first;

	  //
	  // iterate over all points again
	  //

	  for (int j = 0; j < numberPoints; ++j) {
	  
	    //
	    // get j-point handle
	    //

	    const Point & jPoint = _pointData[j].first;

	    //
	    // obtain the correlation Matrix from the
	    // DerivativeCorrelationModel
	    //
	  
	    const Matrix ijCorrelationMatrix = 
	      _correlationModel->getValue(iPoint,
					  jPoint);

	    //
	    // firewall
	    //
	  
	    assert(ijCorrelationMatrix.nrows() == valueDimension);
	    assert(ijCorrelationMatrix.ncols() == valueDimension);

	    //
	    // insert the ijCorrelationMatrix entries into V
	    //

	    for (int k = 0; k < valueDimension; ++k)
	      for (int l = 0; l < valueDimension; ++l)
		V[i + k*numberPoints][j + l*numberPoints] = 
		  ijCorrelationMatrix[k][l];

	  }
	
	}

	//
	// regularize the correlation matrix to alleviate possible poor
	// conditioning
	//

	//       for (int i = 0; i < valueDimension*numberPoints; ++i)
	// 	for (int j = 0; j < valueDimension*numberPoints; ++j)
	// 	  V[i][j] += (10.0 + valueDimension*numberPoints)*
	// 	    std::numeric_limits<double>::epsilon();
      
	return V;

      }

      //
      // compute the value of CorrelationModel at a point
      //

      Matrix
      computeCorrelation(const std::vector<std::pair<Point, Value> > & _pointData,
			 const Point                                 & point,
			 const DerivativeCorrelationModelPointer     & _correlationModel)
      {

	//
	// firewalls
	//

	assert(_pointData.empty() == false);
      
	//
	// get sizes
	//

	const int valueDimension = _pointData.front().second.size();
	const int numberPoints   = _pointData.size();

	//
	// instantiate container
	//
      
	Matrix r(valueDimension*numberPoints,
		 valueDimension);

	//
	// iterate over points
	//

	std::vector<std::pair<Point, Value> >::const_iterator pointDataIter;
	std::vector<std::pair<Point, Value> >::const_iterator pointDataEnd = 
	  _pointData.end();
	int currentPoint = 0;

	for (pointDataIter  = _pointData.begin();
	     pointDataIter != pointDataEnd;
	     ++pointDataIter, ++currentPoint) {
	
	  //
	  // get handle to point
	  //

	  const Point & modelPoint = (*pointDataIter).first;
	
	  //
	  // obtain correlation Matrix from DerivativeCorrelationModel
	  //

	  const Matrix correlationMatrix = 
	    _correlationModel->getValue(modelPoint,
					point);

	  //
	  // insert data into r
	  //

	  for (int i = 0; i < valueDimension; ++i)
	    for (int j = 0; j < valueDimension; ++j)
	      r[currentPoint + i*numberPoints][j] = correlationMatrix[i][j];

	}

	return  r;

      }


      //
      // get values of RegressionModel at all points
      //
    
      Matrix 
      getRegressionModelValues(const std::vector<std::pair<Point, Value> > & _pointData,
			       const DerivativeRegressionModelPointer & _regressionModel)
      {

	//
	// firewalls
	//

	assert(_pointData.empty() == false);

	//
	// get dimensions
	//
      
	const int pointDimension = _pointData.front().first.size();
	const int valueDimension = _pointData.front().second.size();
	const int numberPoints = _pointData.size();
      
	const Dimension regressionModelDimension = 
	  _regressionModel->getDimension(_pointData.front().first);

	assert(regressionModelDimension.first() == 
	       regressionModelDimension.second());

	const int regressionModelDimensionInt = regressionModelDimension.first();

	//
	// firewalls
	//

	assert(pointDimension + 1 == valueDimension);

	//
	// instantiate matrix X
	//

	Matrix X(numberPoints*regressionModelDimensionInt,
		 regressionModelDimensionInt);

	//
	// iterate over points
	//

	std::vector<std::pair<Point, Value> >::const_iterator pointDataIter;
	std::vector<std::pair<Point, Value> >::const_iterator pointDataEnd = 
	  _pointData.end();
	int currentPoint = 0;
      
	for (pointDataIter  = _pointData.begin();
	     pointDataIter != pointDataEnd;
	     ++pointDataIter, ++currentPoint) {
	
	  //
	  // get handle to Point
	  //

	  const Point & point = (*pointDataIter).first;

	  //
	  // get regression value corresponding to point
	  //
	
	  const Matrix regressionModelValues = 
	    _regressionModel->getValues(point);

	  //
	  // copy regressionModelValues into X
	  //

	  for (int i = 0; i < regressionModelDimensionInt; ++i)
	    for (int j = 0; j < regressionModelDimensionInt; ++j)
	      X[currentPoint + i*numberPoints][j] = 
		regressionModelValues[i][j];

	}

	return X;

      }

      //
      // compute matrix A =
      // Inverse[Transpose[X].VInverse.X].Transpose[X].VInverse
      //

      Matrix
      computeMatrixA(Matrix       & _matrixInverseXVX,
		     const Matrix & _matrixX,
		     const Matrix & _matrixInverseV)
      {

	//
	// firewalls
	//

	assert(_matrixX.nrows() == _matrixInverseV.nrows()); // X^T VInverse
	assert(_matrixInverseV.ncols() == _matrixX.nrows());

	//
	// compute X^T
	//

	const Matrix transposeX = transpose(_matrixX);

	//
	// compute X^T VInverse
	//

	const Matrix tranposeXinverseV = mult(transposeX,
					      _matrixInverseV);

	//
	// compute X^T VInverse X
	//

	const Matrix tranposeXinverseVX = mult(tranposeXinverseV,
					       _matrixX);

	//
	// compute X^T VInverse X inverse
	//

	std::pair<Matrix, bool> inverseData = inverse(tranposeXinverseVX);

	assert(inverseData.second == true);
	_matrixInverseXVX = inverseData.first;

	//
	// compute (X^T VInverse X inverse)^(-1) X^T VInverse
	//
      
	return mult(mult(_matrixInverseXVX,
			 transposeX),
		    _matrixInverseV);

      }

      //
      // compute matrix B
      //

      Matrix 
      computeMatrixB(const Matrix & _matrixInverseV,
		     const Matrix & _matrixX,
		     const Matrix & A)
      {

	//
	// firewalls
	//

	assert(_matrixX.ncols() == A.nrows());
	assert(_matrixInverseV.ncols() == _matrixX.nrows());

	//
	// compute X A
	//

	const Matrix XA = mult(_matrixX,
			       A);
      
	assert(XA.nrows() == XA.ncols());

	//
	//
	//

	return mult(_matrixInverseV,
		    identity(XA.nrows()) - XA);

      }

      //
      // compute variance
      //

      double 
      computeVariance(const Vector & _AZ,
		      const Vector & Z,
		      const Matrix & _matrixX,
		      const Matrix & _matrixInverseV,
		      int            numberPoints)
      {

	//
	// compute Z - X.A.Z
	//

	const Vector diffZ = Z - mult(_matrixX,
				      _AZ);

	//
	// compute VInverse.(Z - X.A.Z)
	//

	const Vector vInverseDiffZ = mult(_matrixInverseV,
					  diffZ);

	//
	// compute (scaled) dot product of diffZ and vInverseDiffZ
	//
      
	return mtl::dot(diffZ, vInverseDiffZ)/numberPoints;

      }

      //
      // get max eigenvalue and coresponding eigenvector
      //

      std::pair<double, Vector>
      getSecondPointMomentMaxEigen(const SecondMoment & secondPointMoment)
      {

	//
	// get eigensystem
	//

	const std::pair<Vector, Matrix>  eigenSystem = 
	  secondPointMoment.computeEigenSystem();

	//
	// cache eigenvalues/eigenvectors handles
	//

	const Vector & eigenValues  = eigenSystem.first;
	const Matrix & eigenVectors = eigenSystem.second;

	//
	// find max eigenvalue
	//

	const Vector::const_iterator maxEigenValueIter = 
	  std::max_element(eigenValues.begin(),
			   eigenValues.end());

	//
	// get the max eigenvalue index
	//
      
	const Vector::difference_type maxEigenValueIndex = 
	  std::distance(eigenValues.begin(),
			maxEigenValueIter);
      
	//
	// get a row corresponding to max eigenvalue
	//

	const Matrix::Row eigenVectorRow = 
	  mtl::rows(eigenVectors)[maxEigenValueIndex];


	//
	// copy Matrix::Row into a Vector
	//

	Vector eigenVector;
      
	std::copy(eigenVectorRow.begin(),
		  eigenVectorRow.end(),
		  std::back_insert_iterator<Vector>(eigenVector));
      
	//
	//
	//

	return std::make_pair(*maxEigenValueIter,
			      eigenVector);

      }

      //
      // evaluate the value at point P of the hyper-plane equation with normal
      // n passing through point P0
      //
    
      double 
      evaluateHyperPlaneEquation(const Point  & point,
				 const Vector & normal,
				 const Point  & point0)
      {
      
	return mtl::dot(normal, point - point0);
    
      }
    
    }
  
    //
    // construction/destruction
    //

    DerivativeKrigingModel::DerivativeKrigingModel(const DerivativeRegressionModelPointer & regressionModel,
						   const DerivativeCorrelationModelPointer & correlationModel)
      : _isValid(false),
	_regressionModel(regressionModel),
	_correlationModel(correlationModel)
    {

      return;

    }

    DerivativeKrigingModel::~DerivativeKrigingModel()
    {
    
      return;

    }

    //
    // add point/value pair to the model
    //

    void
    DerivativeKrigingModel::addPoint(const Point & point,
				     const Value & value)
    {

      if (_isValid == false)
	//
	// called when the model is seeded with a collection of points
	//
	_pointData.push_back(std::make_pair(point, value));
      else {

	//
	// this is called when a model already exists and a new point is
	// to be inserted
	//

	_pointData.push_back(std::make_pair(point, value));

	//
	// rebuild model
	//

	build();

      }
      return;

    }

    //
    // get number of points
    //

    int
    DerivativeKrigingModel::getNumberPoints() const
    {

      return _pointData.size();

    }

    //
    // get all points in the kriging model
    //

    std::vector<krigalg::Point> 
    DerivativeKrigingModel::getPoints() const
    {
    
      //
      // instantiate Point container
      //

      std::vector<krigalg::Point> points;
      points.reserve(_pointData.size());
    
      //
      // iterate over _pointData and copy points
      //
    
      std::vector<std::pair<Point, Value> >::const_iterator pointDataIter;
      std::vector<std::pair<Point, Value> >::const_iterator pointDataEnd = 
	_pointData.end();
    
      for (pointDataIter  = _pointData.begin();
	   pointDataIter != pointDataEnd;
	   ++pointDataIter)
	points.push_back((*pointDataIter).first);

      return points;

    }

    //
    // build model using accumulated points
    //

    void
    DerivativeKrigingModel::build()
    {

      // std::cout << "Number points: " << getNumberPoints() << std::endl;

      //
      // reload values at points into a matrix
      //

      const Vector Z = reloadValues(_pointData);

      //
      // construct correlation matrix for all points in the model
      //
    
      const Matrix V = createCorrelationMatrix(_pointData,
					       _correlationModel);
      //
      // compute the inverse of V
      //

      std::pair<Matrix, bool> inverseData = inverse(V);

      assert(inverseData.second == true);
      _matrixInverseV = inverseData.first;

      //
      // get the condition number for V 
      //

      //     std::cout << "Condition number: " 
      // 	      << mtl::one_norm(V)*mtl::one_norm(_matrixInverseV) 
      // 	      << std::endl;

      //
      // get values of RegressionModel at all points
      //
    
      _matrixX = getRegressionModelValues(_pointData,
					  _regressionModel);

      //
      // compute matrix A
      //

      const Matrix A = computeMatrixA(_matrixInverseXVX,
				      _matrixX,
				      _matrixInverseV);

      //
      // compute matrix B
      //

      const Matrix B = computeMatrixB(_matrixInverseV,
				      _matrixX,
				      A);

      //
      // compute AZ and BZ
      //

      _AZ = mult(A,
		 Z);

      _BZ = mult(B, 
		 Z);

      //
      // compute variance
      //

      _sigmaSqr = computeVariance(_AZ,
				  Z,
				  _matrixX,
				  _matrixInverseV,
				  //			_pointData.size());
				  _pointData.size()*
				  (_pointData.front().first.size() + 1));
    
      //
      //
      //

      _isValid = true;

      return;

    }

    //
    // interpolate at point
    //

    Value
    DerivativeKrigingModel::interpolate(const Point & point) const
    {

      //
      // firewalls
      //

      assert(_isValid == true);

      //
      // evaluate RegressionModel at point
      //

      const Matrix Xs = _regressionModel->getValues(point);

      //
      // compute CorrelationModel at point
      //

      const Matrix r = computeCorrelation(_pointData,
					  point,
					  _correlationModel);

      //
      // compute interpolant value at point
      //

      const Vector interpolatedValueVector = 
	mult(Xs, _AZ) + mult(transpose(r), _BZ);		

      //
      // 
      //

      return interpolatedValueVector;

    }

    //
    // get (estimated) interpolation error
    //

    Value
    DerivativeKrigingModel::getMeanSquaredError(const Point & point) const
    {

      //
      // firewalls
      //

      assert(_isValid == true);

      //
      // get the size of Value; since we include the value of the
      // function as well as the gradient the size should equal to the
      // space dimension plus 1
      //

      const int valueDimension = point.size() + 1;

      //
      // evaluate RegressionModel at point
      //

      const Matrix Xs = _regressionModel->getValues(point);

      //
      // compute CorrelationModel at point
      //

      const Matrix r = computeCorrelation(_pointData,
					  point,
					  _correlationModel);

      //
      // compute  u = Xs0 - Transpose[r].VInverse.X
      //

      const Matrix u = Xs - mult(transpose(r), 
				 mult(_matrixInverseV,
				      _matrixX));

      //
      // get self-correlation
      //

      const Matrix sigma = _correlationModel->getValue(point,
						       point);

      //
      // compute the error estimate
      //

      const Matrix errorMatrix = 
	sigma + mult(u, mult(_matrixInverseXVX, transpose(u))) -
	mult(transpose(r), mult(_matrixInverseV, r));
    
      assert(errorMatrix.nrows() == valueDimension);
      assert(errorMatrix.ncols() == valueDimension);

      //
      // extract only self-terms
      //

      Value error(valueDimension);
    
      for (int i = 0; i < valueDimension; ++i)
	error[i] = _sigmaSqr*errorMatrix[i][i];

      //
      //
      //

      return error;

    }
  
    //
    // divide the current DerivativeKrigingModel to create two smaller
    // models; the general idea is to use the second mass moment of all
    // the points involved to detect and bisect the set across its
    // 'largest' dimension
    //

    DerivativeKrigingModel 
    DerivativeKrigingModel::divide()
    {

      //
      // reload points into a container
      //

      std::vector<Point> points;
      points.reserve(_pointData.size());
    
      std::vector<std::pair<Point, Value> >::const_iterator pointDataIter;
      std::vector<std::pair<Point, Value> >::const_iterator pointDataEnd = 
	_pointData.end();

      for(pointDataIter  = _pointData.begin();
	  pointDataIter != pointDataEnd;
	  ++pointDataIter)
	points.push_back((*pointDataIter).first);

      //
      // compute the second mass moment of all points involved
      //
    
      const SecondMoment secondPointMoment(points); 
    
      //
      // use second mass moment ot extract the largest eigenvalue and
      // coresponding to it eigenvector
      //
    
      const std::pair<double, Vector> maxEigenValue = 
	getSecondPointMomentMaxEigen(secondPointMoment);

      const Vector & normal = maxEigenValue.second;

      //
      // obtain handle to the center of mass of all points
      //

      const Point & centerMass = secondPointMoment.getCenterMass();    

      //
      // instantiate childKrigingModel
      //

      DerivativeKrigingModel childKrigingModel(_regressionModel,
					       _correlationModel);

      //
      // iterate over all points; the equation of the plane splitting
      // the model is 
      //
      // f(X) = n.(X - X_c),
      // where: X \in configurational space, X_c - center of mass and 
      // n - the eigenvector of the max second moment eigenvalue
      //
    
      std::vector<std::pair<Point, Value> > pointData;

      for(pointDataIter  = _pointData.begin();
	  pointDataIter != pointDataEnd;
	  ++pointDataIter) {

	//
	// get handle to Point
	//

	const Point & point = (*pointDataIter).first;

	//
	// get handle to Value
	//

	const Value & value = (*pointDataIter).second;

	//
	// points for which the plane equation <= 0 stay here; the rest
	// goes
	//

	if (evaluateHyperPlaneEquation(point, 
				       normal,
				       centerMass) <= 0.0)
	  pointData.push_back(std::make_pair(point,
					     value));
	else
	  childKrigingModel.addPoint(point,
				     value);


      }

      //
      // update and rebuild current kriging model
      //

      _pointData = pointData;
      build();

      //
      // build child kriging model
      //

      childKrigingModel.build();

      //
      //
      //

      return childKrigingModel;

    }

    //
    // output operator
    //

    std::ostream &
    operator<<(std::ostream                 & outputStream,
	       const DerivativeKrigingModel & krigingModel)
    {

      //
      //
      //

      outputStream << "Valid: " << krigingModel._isValid << std::endl;

      //
      // point/value pairs
      //

      std::vector<std::pair<Point, Value> >::const_iterator pointDataIter;
      std::vector<std::pair<Point, Value> >::const_iterator pointDataEnd = 
	krigingModel._pointData.end();

      int currentPointId = 0;

      for (pointDataIter  = krigingModel._pointData.begin();
	   pointDataIter != pointDataEnd;
	   ++pointDataIter, ++currentPointId) {

	outputStream << "Point: " << currentPointId << std::endl;
	outputStream << (*pointDataIter).first << std::endl;
	outputStream << "values: " << std::endl;
	outputStream << (*pointDataIter).second << std::endl;

      }

      //
      // arrays
      //

      outputStream << "matrixX: " << std::endl;
      outputStream << krigingModel._matrixX << std::endl;

      outputStream << "matrixInverseV: " << std::endl;
      outputStream << krigingModel._matrixInverseV << std::endl;

      outputStream << "matrixInverseXVX: " << std::endl;
      outputStream << krigingModel._matrixInverseXVX << std::endl;

      outputStream << "AZ: " << std::endl;
      outputStream << krigingModel._AZ << std::endl;

      outputStream << "BZ: " << std::endl;
      outputStream << krigingModel._BZ << std::endl;

      //
      // sigmaSqr
      //

      outputStream << "sigmaSqr: " << krigingModel._sigmaSqr << std::endl;

      //
      // regression and correlation models
      //

      //    outputStream << "regressionModel: " << *(krigingModel._regressionModel)
      // 		 << std::endl;

      //     outputStream << "correlationModel: " << *(krigingModel._correlationModel)
      // 		 << std::endl;

      return outputStream;

    }

}


