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
// File:        MultivariateKrigingModel.cc
// Package:     kriging algorithm
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: Class implementing kriging model including derivative 
// Description: information.
//

#include "MultivariateKrigingModel.h"

#include "CorrelationModelFactory.h"
#include "RegressionModelFactory.h"
#include "SecondMoment.h"

#include "toolbox/database/Database.h"
#include "toolbox/base/Utilities.h"

#include <mtl/mtl.h>
#include <mtl/mtl2lapack.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <limits>

//
//
// 

namespace krigalg {

    //
    //
    //

    namespace {

      //
      // const data
      //

      const std::string numberPointsKey("number_points");
      const std::string numberValuesKey("number_value");
      const std::string pointDimensionKey("point_dimension");
      const std::string pointDataKey("point_data");
      const std::string pointValuesKey("point_values");
      const std::string regressionModelClassNameKey("regression_model_class_name");
      const std::string correlationModelClassNameKey("correlation_model_class_name");

    
      //
      // construct correlation matrix for all points in the model
      //   
    
      Matrix 
      createCorrelationMatrix(const std::vector<Point> & _points,
			      CorrelationModelPointer  & _correlationModel,
			      int                         pointDimension,
			      int                         numberPoints)
      {

	//
	// firewalls
	//

	assert(_points.empty() == false);
	assert(_points.size()  == numberPoints);

	//
	// instantiate matrix
	//

	Matrix V(numberPoints,
		 numberPoints);

	//
	// iterate over all points
	//

	for (int i = 0; i < numberPoints; ++i) {

	  //
	  // get i-point handle 
	  //

	  const Point & iPoint = _points[i];

	  //
	  // iterate over all points again
	  //

	  for (int j = 0; j < numberPoints; ++j) {
	  
	    //
	    // get j-point handle
	    //

	    const Point & jPoint = _points[j];

	    //
	    // obtain the value of the correlation from CorrelationModel
	    //

	    V[i][j] = (_correlationModel->getValue(iPoint,
						   jPoint))[0][0];

	  }
	
	}

	//
	// regularize the correlation matrix to alleviate possible poor
	// conditioning
	//

	for (int i = 0; i < numberPoints; ++i)
	  V[i][i] += (10.0 + numberPoints)*
	    std::numeric_limits<double>::epsilon();
      
	return V;

      }
    
      //
      // get values of RegressionModel at all points
      //
    
      Matrix 
      getRegressionModelValues(const std::vector<Point>     & _points,
			       const RegressionModelPointer & _regressionModel,
			       int                             pointDimension,
			       int                             numberPoints)
      {

	//
	// firewalls
	//

	assert(_points.empty() == false);
	assert(_points.size()  == numberPoints);

	//
	// get dimensions
	//

	const Dimension regressionModelDimensionTemp = 
	  _regressionModel->getDimension(_points.front());

	assert(regressionModelDimensionTemp.second() == 1);
      
	const int regressionModelDimension = 
	  regressionModelDimensionTemp.first();

	//
	// instantiate matrix X
	//

	Matrix X(numberPoints,
		 regressionModelDimension);

	//
	// iterate over points
	//

	std::vector<Point>::const_iterator pointsIter;
	std::vector<Point>::const_iterator pointsEnd = _points.end();
	int currentPoint = 0;
      
	for (pointsIter  = _points.begin();
	     pointsIter != pointsEnd;
	     ++pointsIter, ++currentPoint) {
	
	  //
	  // get handle to Point
	  //

	  const Point & point = *pointsIter;

	  //
	  // get regression value corresponding to point
	  //

	  const Matrix regressionModelValuesTemp = 
	    _regressionModel->getValues(point);
	
	  assert(regressionModelValuesTemp.ncols() == 1);

	  Vector regressionModelValues(regressionModelValuesTemp.nrows());
	  
	  for (int i = 0; i < regressionModelValuesTemp.nrows(); ++i)
	    regressionModelValues[i] = regressionModelValuesTemp[i][0];

	  //
	  // copy regressionModelValues into X
	  //

	  std::copy(regressionModelValues.begin(),
		    regressionModelValues.end(),
		    &(X[currentPoint][0]));

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
	// compute X^T VInverse
	//

	const Matrix tranposeXinverseV = mult(_matrixX,
					      _matrixInverseV,
					      true,
					      false);

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
			 _matrixX,
			 false,
			 true),
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
      // reload values into Z
      //

      std::vector<Vector>
      reloadValues(const std::vector<std::vector<Value> > & _values,
		   int                                      pointDimension,
		   int                                      numberPoints,
		   int                                      numberValues)
      {

	//
	// allocate space for vectors Z
	//
      
	std::vector<Vector> Z;
      
	for (int iValue = 0; iValue < numberValues; ++iValue)
	  Z.push_back(Vector(numberPoints));
      
	std::vector<std::vector<Value> >::const_iterator valuesIter;
	std::vector<std::vector<Value> >::const_iterator valuesEnd = 
	  _values.end();
	int currentPoint = 0;

	for (valuesIter  = _values.begin();
	     valuesIter != valuesEnd;
	     ++valuesIter, ++currentPoint) {
	
	  //
	  // get handle to all values coresponding to a point
	  //
      
	  const std::vector<Value> & pointValues = *valuesIter;
	
	  //
	  // firewalls
	  //
	
	  assert(pointValues.size() == numberValues);
	
	  for (int iValue = 0; iValue < numberValues; ++iValue) {

	    //
	    // get handle to value
	    //

	    const Value & value = pointValues[iValue];

	    //
	    // copy Value
	    //

	    for (Value::size_type i = 0; i < value.size(); ++i)
	      Z[iValue][currentPoint + i*numberPoints] = value[i];

	  }
	
	}

	return Z;

      }

      //
      // compute _AZ and _BZ
      //

      std::pair<std::vector<Vector>, std::vector<Vector> >
      computeAZBZPair(const std::vector<Vector> & Z,
		      const Matrix              & A,
		      const Matrix              & B)
      {
      
	//
	// allocate space for _AZ and _BZ
	//
      
	std::vector<Vector> _AZ;
	std::vector<Vector> _BZ;

	//
	// iterate through Z
	//

	std::vector<Vector>::const_iterator zIter;
	const std::vector<Vector>::const_iterator zEnd = Z.end();

	for (zIter  = Z.begin(); zIter != zEnd; ++zIter) {

	  //
	  // get handle to part of Z coresponding to value
	  //

	  const Vector & zValue = *zIter;

	  //
	  // multiply and store
	  //

	  _AZ.push_back(mult(A, zValue));
	  _BZ.push_back(mult(B, zValue));

	}

	//
	//
	//

	return std::make_pair(_AZ, _BZ);

      }
    
      //
      // compute variance
      //

      std::vector<double >
      computeVariance(const std::vector<Vector> & _AZ,
		      const std::vector<Vector> & Z,
		      const Matrix              & _matrixX,
		      const Matrix              & _matrixInverseV,
		      int                         numberPoints,
		      int                         numberPerPointValues)
      {

	//
	// allocate sigmaSqr
	//

	std::vector<double> sigmaSqr;
	sigmaSqr.reserve(Z.size());

	//
	// iterate over all values
	//

	for (int iValue = 0; iValue < numberPerPointValues; ++iValue) {
	
	  //
	  // compute Z - X.A.Z
	  //
	
	  const Vector diffZ = Z[iValue] - mult(_matrixX,
						_AZ[iValue]);
	
	  //
	  // compute VInverse.(Z - X.A.Z)
	  //
	
	  const Vector vInverseDiffZ = mult(_matrixInverseV,
					    diffZ);
	
	  //
	  // compute (scaled) dot product of diffZ and vInverseDiffZ
	  //
	
	  sigmaSqr.push_back(dot(diffZ, vInverseDiffZ)/numberPoints);

	}

	return sigmaSqr;
	
      }

      //
      // compute the value of CorrelationModel at a point
      //

      Vector
      computeCorrelation(const std::vector<Point>      & _points,
			 const Point                   & point,
			 const CorrelationModelPointer & _correlationModel)
      {

	//
	// firewalls
	//

	assert(_points.empty() == false);
      
	//
	// get sizes
	//

	const int numberPoints   = _points.size();

	//
	// instantiate container
	//
      
	Vector r(numberPoints);

	//
	// iterate over points
	//

	std::vector<Point>::const_iterator pointsIter;
	std::vector<Point>::const_iterator pointsEnd = _points.end();
	int currentPoint = 0;

	for (pointsIter  = _points.begin();
	     pointsIter != pointsEnd;
	     ++pointsIter, ++currentPoint) {
	
	  //
	  // get handle to point
	  //

	  const Point & modelPoint = *pointsIter;
	
	  //
	  // obtain correlation Matrix from CorrelationModel
	  //

	  r[currentPoint] = (_correlationModel->getValue(modelPoint,
							 point))[0][0];

	}

	return  r;

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
      
	return dot(normal, point - point0);
    
      }

      //
      // copy i-th row data
      //

      inline void
      getRow(Vector       & row,
	     const Matrix & matrix,
	     int            i)
      {

	//
	// firewalls
	//

	assert(row.size() >= matrix.ncols());
	assert((int)Matrix::orientation::id == (int)mtl::ROW_MAJOR);

	//
	//
	//

	std::copy(matrix[i].begin(),
		  matrix[i].end(),
		  &(row[0]));

	//
	//
	//

	return;

      }

      //
      // copy i-th column data
      //

      inline void
      getColumn(Vector       & column,
		const Matrix & matrix,
		int            i)
      {
	
	//
	// firewalls
	//

	assert(column.size() >= matrix.nrows());

	//
	//
	//
	
	for (int j = 0; j < matrix.nrows(); ++j)
	  column[j] = matrix[j][i];


	//
	//
	//

	return;

      }

      //
      // pack a contents of a matrix into a vector container
      //

      inline void
      packMatrix(std::vector<double> & packedContainer,
		 const Matrix        & matrix)
      {

	const int numberColumns = matrix.ncols();
	const int numberRows    = matrix.nrows(); 
	
	const double * data = matrix.data();
	
	//
	// store order
	//

	packedContainer.push_back(numberRows);
	packedContainer.push_back(numberColumns);

	//
	//
	//

	packedContainer.insert(packedContainer.end(),
			       &(data[0]),
			       &(data[numberRows*numberColumns]));
	
	//
	//
	//

	return;

      }

      //
      // unpack Matrix
      //

      inline Matrix
      unpackMatrix(const std::vector<double> & packedContainer,
		   int                       & currentOffset)
      {

	//
	// get order
	//

	const int numberRows    = (int)packedContainer[currentOffset++];
	const int numberColumns = (int)packedContainer[currentOffset++];

	//
	// instantiate matrix
	//

	Matrix matrix(numberRows,
		      numberColumns);
	
	//
	// fill in 
	//

	const int numberData = numberRows*numberColumns;

	std::copy(&(packedContainer[currentOffset]),
		  &(packedContainer[currentOffset + numberData]),
		  matrix.data());
	
	currentOffset += numberData;

	//
	//
	//
	
	return matrix;

      }
		   
            
      //
      // pack std::vector<Vector>
      //

      inline void
      packVectorVector(std::vector<double>       & packedContainer,
		       const std::vector<Vector> & vectorVector)
      {

	std::vector<Vector>::const_iterator vectorVectorIter;
	std::vector<Vector>::const_iterator vectorVectorEnd = 
	  vectorVector.end();
	
	packedContainer.push_back(vectorVector.size());

	for (vectorVectorIter  = vectorVector.begin();
	     vectorVectorIter != vectorVectorEnd;
	     ++vectorVectorIter) {

	  //
	  // store size
	  //
	  
	  packedContainer.push_back((*vectorVectorIter).size());

	  //
	  // store data
	  //

	  packedContainer.insert(packedContainer.end(),
				 (*vectorVectorIter).begin(),
				 (*vectorVectorIter).end());
	}

	//
	//
	//

	return;

      }

      //
      // unpack std::vector<Vector>
      //

      inline std::vector<Vector>
      unpackVectorVector(const std::vector<double> & packedContainer,
			 int                       & currentOffset)
      {

	//
	// extract std::vector<Vector> size
	//

	const int vectorVectorSize = (int)packedContainer[currentOffset++];

	//
	// instantiate container
	//

	std::vector<Vector> vectorVector;
	vectorVector.reserve(vectorVectorSize);

	//
	// iterate over all std::vector<Vector> entries
	//

	for (int i = 0; i < vectorVectorSize; ++i) {

	  //
	  // get size
	  //
	  
	  const int vectorSize = (int)packedContainer[currentOffset++];

	  //
	  // instantiate Vector
	  //

	  const Vector vector(vectorSize,
			      &(packedContainer[currentOffset]));

	  
	  currentOffset += vectorSize;

	  //
	  // insert vector into vectorVector
	  //

	  vectorVector.push_back(vector);

	}

	return vectorVector;

      }
      
      //
      // Check if for the current collection of points their
      // correlation matrix is well conditioned. Well conditioned
      // matrix is indicated by true and poorly conditioned matrix by
      // false.
      //

      inline bool
      checkCorrelationMatrixConditionNumber(const Point              & point,
					    const std::vector<Point> & _points,
					    CorrelationModelPointer  & _correlationModel,
					    int pointDimension)
      {

	//
	// construct correlation matrix for all points in the model
	//
    
	std::vector<Point> pointsCopy(_points.begin(),
				      _points.end());
	pointsCopy.push_back(point);

	const Matrix V = createCorrelationMatrix(pointsCopy,
						 _correlationModel,
						 pointDimension,
						 pointsCopy.size());

	//
	// compute the inverse of V
	//
    
	std::pair<Matrix, bool> inverseData = inverse(V);

	//
	// should inversion fail return false
	//
	if (inverseData.second == false)
	  return false;

	const Matrix matrixInverseV = inverseData.first;

	//
	// get the condition number for V 
	//

	const double vConditionNumber = 
	  mtl::one_norm(V)*mtl::one_norm(matrixInverseV);

	const double maxConditionNumber = 1.0e10;
	
	if (vConditionNumber > maxConditionNumber)
	  return false;

	
	return true;

      }

      Vector
      getPointExtent(const std::vector<Point> & points)
      {

        //
        // firewalls
        //

        assert(points.empty() == false);

        //
        //
        //

        const int pointDimension = points.front().size();

        //
        // instantiate containers to store the min/max
        //

        Vector minExtent(pointDimension,
                         std::numeric_limits<double>::max());

        Vector maxExtent(pointDimension,
                         std::numeric_limits<double>::min());

        //
        // iterate over points
        //

        std::vector<Point>::const_iterator pointsIter;
        const std::vector<Point>::const_iterator pointsEnd = points.end();

        for (pointsIter  = points.begin();
             pointsIter != pointsEnd;
             ++pointsIter) {

          //
          //
          //

          const Vector & point = *pointsIter;

          //
          // iterate over dimensions
	  //
	  
          for (int i = 0; i < pointDimension; ++i) {
            minExtent[i] = std::min(minExtent[i], point[i]);
            maxExtent[i] = std::max(maxExtent[i], point[i]);
          }

        }

        return maxExtent - minExtent;

      }

      //
      // given a Vector compute maximum value of all its components
      //
      
      inline double
      getMaxVectorComponent(const Vector & vector)
      {

        return *(std::max_element(vector.begin(),
                                  vector.end()));

      }

    }

    //
    // construction/destruction
    //
    MultivariateKrigingModel::MultivariateKrigingModel(const RegressionModelPointer & regressionModel,
						       const CorrelationModelPointer & correlationModel)
      :  _isValid(false),
	 _regressionModel(regressionModel),
	 _correlationModel(correlationModel)
    {

      return;

    }

    MultivariateKrigingModel::~MultivariateKrigingModel()
    {

      return;

    }

    //
    //
    //

    bool
    MultivariateKrigingModel::hasError() const
    {
      return true;
    }

    //
    // 
    //

    bool
    MultivariateKrigingModel::hasGradient() const
    {
      return false;
    }
    
    //
    // clone
    //

    MultivariateKrigingModel *
    MultivariateKrigingModel::clone() const
    {

      return new MultivariateKrigingModel(*this);

    }


    //
    // add point and values to the model
    //
  
    bool 
    MultivariateKrigingModel::addPoint(const Point              & point,
				       const std::vector<Value> & values)
    {

      if (_isValid == false) {

	//
	// firewalls
	//

	assert(values.empty() == false);
	assert(values.front().size() == 1);

	//
	// check if insertion of a point would cause conditioning problems
	//

	if (checkCorrelationMatrixConditionNumber(point,
						  _points,
						  _correlationModel,
						  getPointDimension()) == false)
	  return false;
	
	//
	// insert into containers
	//
	
	_points.push_back(point);
	_values.push_back(values);

	//
	// attempt to build model
	//
	
	build();

      } else {
      
	//
	// this gets executed when the model already exists and a new point is
	// to be inserted
	//
	
	if (checkCorrelationMatrixConditionNumber(point,
						  _points,
						  _correlationModel,
						  getPointDimension()) == false)
	  return false;


	//
	// insert into containers
	//
	
	_points.push_back(point);
	_values.push_back(values);
	
	//
	// re-build the model
	//

	build();

      }

      //
      //
      //

      return true;

    }

    //
    // get number of points in the model
    //
  
    int
    MultivariateKrigingModel::getNumberPoints() const
    {

      return _points.size();

    }
  
    //
    // get handle to kriging model points
    //

    const std::vector<Point> &
    MultivariateKrigingModel::getPoints() const
    {

      return _points;

    }

    //
    // get number of values stored at a point
    //
  
    int 
    MultivariateKrigingModel::getNumberValues() const
    {

      if (_values.empty() == true)
	return 0;
      else 
	return _values.front().size();

    }

    //
    // get point dimension
    //
  
    int 
    MultivariateKrigingModel::getPointDimension() const
    {

      if (_points.empty() == true) 
	return 0;
      else
	return _points.front().size();

    }
  
    //
    // get value dimension
    //
  
    int 
    MultivariateKrigingModel::getValueDimension() const
    {

      if (_values.empty() == true)
	return 0;
      else 
	return _values.front().front().size();

    }

    //
    // get regression model
    //

    RegressionModelPointer
    MultivariateKrigingModel::getRegressionModel() const
    {

      return _regressionModel;

    }

    //
    // get correlation model
    //

    CorrelationModelPointer
    MultivariateKrigingModel::getCorrelationModel() const
    {

      return _correlationModel;
    
    }


    //
    // build the model using accumulated points
    //
  
    void 
    MultivariateKrigingModel::build()
    {
    
      //
      // get dimensions
      //

      const int pointDimension = getPointDimension();
      const int numberPoints   = getNumberPoints();
      const int numberValues   = getNumberValues();

      //
      // do not build the model unless there are at least min number
      // of points
      //
      std::cout << "Number points: " << numberPoints << std::endl;
      if (numberPoints <= pointDimension)
	return;
      std::cout << "Building a full model" << std::endl;
      //
      // firewall
      //

      assert(_points.size() == _values.size());

      //
      // generate value-independent data
      //
#if 0      
      const double dMax = getMaxVectorComponent(getPointExtent(_points))*5;
      const double theta = 0.5/(dMax*dMax);
      std::cout << "setting theta to " << theta << std::endl;
      _correlationModel->setThetas(std::vector<double>(1, theta));
#endif 
      //
      // construct correlation matrix for all points in the model
      //
    
      const Matrix V = createCorrelationMatrix(_points,
					       _correlationModel,
					       pointDimension,
					       numberPoints);

      //
      // compute the inverse of V
      //
    
      std::pair<Matrix, bool> inverseData = inverse(V);

      assert(inverseData.second == true);
      _matrixInverseV = inverseData.first;

      //
      // get the condition number for V 
      //

//       std::cout << "Condition number: " 
// 		<< mtl::one_norm(V)*mtl::one_norm(_matrixInverseV) 
// 		<< std::endl;
    
      //
      // get values of RegressionModel at all points
      //
    
      _matrixX = getRegressionModelValues(_points,
					  _regressionModel,
					  pointDimension,
					  numberPoints);

      //
      // compute V^{-1} X 
      //

      _matrixInverseVX = mult(_matrixInverseV,
			      _matrixX);

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
      // compute all value dependent data
      //
    
      //
      // compute reloaded data
      //

      const std::vector<Vector> Z = reloadValues(_values,
						 pointDimension,
						 numberPoints,
						 numberValues);
    
      //
      // compute _AZ and _BZ
      //

      std::pair<std::vector<Vector>, std::vector<Vector> > aZbZPair = 
	computeAZBZPair(Z,
			A,
			B);

      _AZ = aZbZPair.first;
      _BZ = aZbZPair.second;

      //
      // compute variance
      //

      _sigmaSqr = computeVariance(_AZ,
				  Z,
				  _matrixX,
				  _matrixInverseV,
				  numberPoints,
				  numberValues);
    
      //
      // make model valid
      //

      _isValid = true;

      //
      //
      // 

      return;

    }

    //
    // check if model is valid
    //

    bool
    MultivariateKrigingModel::isValid() const
    {

      return _isValid;

    }

    //
    // interpolate value at a point
    //

    Value
    MultivariateKrigingModel::interpolate(int           valueId,
					  const Point & point) const
    {

      //
      // firewalls
      //

      assert(_isValid == true);

      //
      // evaluate RegressionModel at point
      //

      const Matrix XsTemp = _regressionModel->getValues(point);

      assert(XsTemp.ncols() == 1);

      Vector Xs(XsTemp.nrows());

      for (int i = 0; i < XsTemp.nrows(); ++i)
	Xs[i] = XsTemp[i][0];

      //
      // compute CorrelationModel at point
      //

      const Vector r = computeCorrelation(_points,
					  point,
					  _correlationModel);

      //
      // compute interpolant value at point
      //

      const Vector interpolatedValueVector(1, 
					   dot(Xs, _AZ[valueId]) + dot(r, _BZ[valueId]));

      //
      //
      //

      assert(interpolatedValueVector.size() == 1);

      //
      // 
      //

      return interpolatedValueVector;

    }

    //
    // get (estimated) interpolation error
    //

    Value
    MultivariateKrigingModel::getMeanSquaredError(int           valueId,
						  const Point & point) const
    {

      //
      // firewalls
      //

      assert(_isValid == true);

      //
      // evaluate RegressionModel at point
      //

      const Matrix XsTemp = _regressionModel->getValues(point);

      assert(XsTemp.ncols() == 1);

      Vector Xs(XsTemp.nrows());

      for (int i = 0; i < XsTemp.nrows(); ++i)
	Xs[i] = XsTemp[i][0];

      //
      // compute CorrelationModel at point
      //

      const Vector r = computeCorrelation(_points,
					  point,
					  _correlationModel);
      
      //
      // compute  u = Xs0 - Transpose[r].VInverse.X
      //

      const Vector u = Xs - mult(transpose(_matrixInverseVX),
				 r);

      //
      // get self-correlation
      //

      const double sigma = (_correlationModel->getValue(point,
							point))[0][0];

      //
      // compute the error estimate
      //

      Vector errorVector(1);
      double errorValue;

      //
      // initialize with the self correlation term
      //
      
      errorValue = sigma;

      //
      // add u.(XVX)^-1.u^T contribution
      //
      
      errorValue += dot(u, 
			mult(_matrixInverseXVX, u));
	
      //
      // add r^T V^-1 r contribution
      //
      
      errorValue -= dot(r,
			mult(_matrixInverseV, r));
      
      //
      // apply self-correlation
      //
      
      errorValue *= _sigmaSqr[valueId];

      //
      //
      //

      errorVector[0] = errorValue;

      //
      //
      //

      return errorVector;

    }  

    //
    // divide the current KrigingModel to create two smaller
    // models; the general idea is to use the second mass moment of all
    // the points involved to detect and bisect the set across its
    // 'largest' dimension
    //

    MultivariateKrigingModel 
    MultivariateKrigingModel::divide()
    {

      //
      // compute the second mass moment of all points involved
      //
    
      const SecondMoment secondPointMoment(_points); 
    
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

      MultivariateKrigingModel childKrigingModel(_regressionModel,
							   _correlationModel);

      //
      // iterate over all points; the equation of the plane splitting
      // the model is 
      //
      // f(X) = n.(X - X_c),
      // where: X \in configurational space, X_c - center of mass and 
      // n - the eigenvector of the max second moment eigenvalue
      //
    
      std::vector<Point>               points;
      std::vector<std::vector<Value> > values;

      std::vector<Point>::size_type iPoint;
      std::vector<Point>::size_type numberPoints = _points.size();

      for (iPoint = 0; iPoint < numberPoints; ++iPoint) {

	//
	// get handle to Point
	//

	const Point & point = _points[iPoint];

	//
	// get handle to Value
	//

	const std::vector<Value> & value = _values[iPoint];

	//
	// points for which the plane equation <= 0 stay here; the rest
	// goes
	//

	if (evaluateHyperPlaneEquation(point, 
				       normal,
				       centerMass) <= 0.0) {

	  points.push_back(point);
	  values.push_back(value);
	
	} else
	  childKrigingModel.addPoint(point,
				     value);


      }

      //
      // update and rebuild current kriging model (make sure that 
      // it makes sense for the model to be rebuilt)
      //

      _points = points;
      _values = values;
      build();

      //
      // build child kriging model (only if this makes sense)
      //

      childKrigingModel.build();

      //
      //
      //

      return childKrigingModel;

    }

    //
    // output object data to database; store only point-value pairs and
    // rely in the getFromDatabase() to properly build the model
    //

    void
    MultivariateKrigingModel::putToDatabase(toolbox::Database& db) const
    {

      //
      // firewalls
      //

      assert(_isValid == true);
      assert(_points.size() == _values.size());

      //
      // 
      // 

      const int numberPoints        = getNumberPoints();
      const int numberPointValues   = getNumberValues();
      const int pointDimension      = getPointDimension();

      //
      // number of point/value pairs
      //

      db.putInteger(numberPointsKey, numberPoints);

      //
      // number of values at a point
      //

      db.putInteger(numberValuesKey, numberPointValues);

      //
      // store point dimension
      //

      db.putInteger(pointDimensionKey, pointDimension);

      //
      // iterate through point/value pairs
      //

      std::vector<double> allPointData;
      std::vector<double> allValueData;
      

      for(std::vector<Point>::size_type i = 0; i < numberPoints; ++i) {

	//
	// store point
	//

	const Point & point = _points[i];

	allPointData.insert(allPointData.end(),
			    point.begin(),
			    point.end());

	//
	// store values
	//
	
	const std::vector<Value> & pointValues = _values[i];
	
	for (std::vector<Value>::size_type j = 0; j < numberPointValues; ++j)
	  allValueData.insert(allValueData.end(),
			      pointValues[j].begin(),
			      pointValues[j].end());

      }

      db.putDoubleArray(pointDataKey, allPointData);
      db.putDoubleArray(pointValuesKey, allValueData);

      //
      // store regression model
      //

      db.putString(regressionModelClassNameKey, 
		   _regressionModel->getClassName());

      _regressionModel->putToDatabase(db);

      //
      // store correlation model
      //

      db.putString(correlationModelClassNameKey,
		   _correlationModel->getClassName());

      _correlationModel->putToDatabase(db);
	
      //
      //
      //

      return;

    }

    //
    // fill in the object from Database
    //

    void
    MultivariateKrigingModel::getFromDatabase(toolbox::Database& db)
    {

      //
      // get number point/value pairs
      //

      const int numberPoints = db.getInteger(numberPointsKey);

      //
      // get number vales at a point
      //

      const int numberValues = db.getInteger(numberValuesKey);

      //
      // get point dimension
      //
      
      const int pointDimension = db.getInteger(pointDimensionKey);

      //
      // get point/value data
      //

      std::vector<double> allPointData;
      db.getDoubleArray(pointDataKey, allPointData);
      
      std::vector<double> allValueData;
      db.getDoubleArray(pointValuesKey, allValueData);

      //
      // firewalls
      //

      assert(allPointData.size() == numberPoints*pointDimension);
      assert(allValueData.size() == numberPoints*numberValues);
      
      //
      // reload point/value data
      //

      for (int i = 0; i < numberPoints; ++i) {

	//
	// instantiate and store point data 
	//
	
	Point point(pointDimension,
		    &(allPointData[i*pointDimension]));
	
	_points.push_back(point);
	
	//
	// get values coresponding to point data
	//
		
	std::vector<Value> pointValues;

	for (int j = 0; j < numberValues; ++j) {

	  //
	  // instantiate and create value
	  //

	  Value value(1, 
		      &(allValueData[i*numberValues +
				     j]));

	  //
	  // store value
	  //
	  
	  pointValues.push_back(value);

	}
    
	//
	// store point values
	//

	_values.push_back(pointValues);

      }

      //
      // get the class name of the regression model
      //

      const std::string regressionModelClassName = 
	db.getString(regressionModelClassNameKey);

      //
      // build the regression model
      //

      _regressionModel =
	RegressionModelFactory().build(regressionModelClassName);

      //
      // load regression model data from db
      //

      _regressionModel->getFromDatabase(db);

      //
      // get the class name for the correlation model
      //

      const std::string correlationModelClassName = 
	db.getString(correlationModelClassNameKey);

      //
      // build correlation model
      //

      _correlationModel =
	CorrelationModelFactory().build(correlationModelClassName);

      //
      // load correlation model data from db
      //

      _correlationModel->getFromDatabase(db);

      //
      // rebuild model
      //
      
      build();

      //
      //
      //

      return;

    }

    //
    // pack data into an array of doubles
    //

    void
    MultivariateKrigingModel::pack(std::vector<double> & packedContainer) const 
    {

      //
      // firewalls
      //
      
      assert(_isValid == true);
      assert(_points.size() == _values.size());

      //
      // clear container
      //

      packedContainer.clear();

      //
      // store size data
      //

      const int numberPoints      = getNumberPoints();
      const int numberPointValues = getNumberValues();

      packedContainer.push_back(numberPoints);
      packedContainer.push_back(numberPointValues);
      packedContainer.push_back(getPointDimension());
      packedContainer.push_back(getValueDimension());

      //
      // store point/value data
      //

      for(std::vector<Point>::size_type i = 0; i < numberPoints; ++i) {
	
	//
	// store point
	//

	const Point & point = _points[i];

	packedContainer.insert(packedContainer.end(),
			       point.begin(),
			       point.end());

	//
	// store values
	//
	
	const std::vector<Value> & pointValues = _values[i];
	
	for (std::vector<Value>::size_type j = 0; j < numberPointValues; ++j)
	  packedContainer.insert(packedContainer.end(),
				 pointValues[j].begin(),
				 pointValues[j].end());
	
      }

      //
      // _matrixX
      //

      packMatrix(packedContainer,
		 _matrixX);
      
      //
      // _matrixInverseV
      //
      
      packMatrix(packedContainer,
		 _matrixInverseV);

      //
      // _matrixInverseXVX
      //

      packMatrix(packedContainer,
		 _matrixInverseXVX);
      
      //
      // _matrixInverseVX
      //

      packMatrix(packedContainer,
		 _matrixInverseVX);

      //
      // _AZ
      //

      packVectorVector(packedContainer,
		       _AZ);
      
      //
      // _BZ
      //

      packVectorVector(packedContainer,
		       _BZ);

      //
      // _sigmaSqr
      //

      packedContainer.push_back(_sigmaSqr.size());
      
      packedContainer.insert(packedContainer.end(),
			     _sigmaSqr.begin(),
			     _sigmaSqr.end());

      //
      // regression model
      // 

      RegressionModelFactory::ClassIdentifier regressionClassId = 
	RegressionModelFactory::getClassId(_regressionModel->getClassName());
      packedContainer.push_back(regressionClassId);

      //
      // correlation model
      //

      CorrelationModelFactory::ClassIdentifier correlationClassId = 
	CorrelationModelFactory::getClassId(_correlationModel->getClassName());

      packedContainer.push_back(correlationClassId);

      std::vector<double> thetas;
      _correlationModel->getThetas(thetas);

      packedContainer.push_back(thetas.size());
      packedContainer.insert(packedContainer.end(),
			     thetas.begin(),
			     thetas.end());

      //
      //
      //

      return;

    }

    //
    // initialize data from an array of double
    //

    void
    MultivariateKrigingModel::unpack(const std::vector<double> & packedContainer)
    {

      //
      // 
      //

      int currentOffset = 0;
      
      //
      // get size data
      //

      const int numberPoints      = (int)packedContainer[currentOffset++];
      const int numberPointValues = (int)packedContainer[currentOffset++];
      const int pointDimension    = (int)packedContainer[currentOffset++];

      //
      // clear _points and _values
      //

      _points.clear();
      _values.clear();

      _points.reserve(numberPoints);
      _points.reserve(numberPoints);

      //
      // get point value data
      //

      for (int i = 0; i < numberPoints; ++i) {

	//
	// instantiate point
	//

	const Point point(pointDimension,
			  &(packedContainer[currentOffset]));

	currentOffset += pointDimension;

	//
	// insert point
	//

	_points.push_back(point);

	//
	// instantiate values
	//

	std::vector<Value> pointValues;
	pointValues.reserve(numberPointValues);

	for (int j = 0; j < numberPointValues; ++j) {

	  //
	  // instantiate value
	  //

	  Vector value(1,
		       &(packedContainer[currentOffset]));
	  
	  ++currentOffset;
	  
	  //
	  // insert value
	  //
	  
	  pointValues.push_back(value);

	}

	//
	// store pointValues
	//

	_values.push_back(pointValues);

      }
	
      //
      // unpack _matrixX
      //
      
      _matrixX = unpackMatrix(packedContainer,
			      currentOffset);

      //
      // unpack _matrixInverseV
      //

      _matrixInverseV = unpackMatrix(packedContainer,
				     currentOffset);

      //
      // unpack _matrixInverseXVX
      //

      _matrixInverseXVX = unpackMatrix(packedContainer,
				       currentOffset);

      //
      // unpack _matrixInverseVX
      //

      _matrixInverseVX = unpackMatrix(packedContainer,
				      currentOffset);

      //
      // unpack _AZ
      //

      _AZ = unpackVectorVector(packedContainer,
			       currentOffset);

      //
      // unpack _BZ
      //

      _BZ = unpackVectorVector(packedContainer,
			       currentOffset);

      //
      // unpack _sigmaSqr
      //

      const int sigmaSqrSize = (int)packedContainer[currentOffset++];

      _sigmaSqr = 
	std::vector<double>(&(packedContainer[currentOffset]),
			    &(packedContainer[currentOffset + sigmaSqrSize]));

      currentOffset += sigmaSqrSize;

      //
      // unpack regression model
      //

      const RegressionModelFactory::ClassIdentifier 
	regressionClassId = (RegressionModelFactory::ClassIdentifier)packedContainer[currentOffset++];

      _regressionModel = 
	RegressionModelFactory().build(regressionClassId);

      //
      // unpack correlation model
      //

      const CorrelationModelFactory::ClassIdentifier
	correlationClassId = (CorrelationModelFactory::ClassIdentifier)packedContainer[currentOffset++];
      const int thetaSize = (int)packedContainer[currentOffset++];
      
      const std::vector<double> thetas(&(packedContainer[currentOffset]),
				       &(packedContainer[currentOffset +
							 thetaSize]));

      currentOffset += thetaSize;
      
      _correlationModel = 
	CorrelationModelFactory().build(correlationClassId,
						  thetas);


      //
      // make valid
      //

      _isValid = true;

      //
      //
      //

      return;

    }

    //
    // output operator
    //

    std::ostream &
    operator<<(std::ostream                   & outputStream,
	       const MultivariateKrigingModel & krigingModel)
    {

      //
      //
      //

      outputStream << "Valid: " << krigingModel._isValid << std::endl;

      //
      // point/value pairs
      //


      std::vector<Point>::size_type iPoint;
      const std::vector<Point>::size_type numberPoints = 
	krigingModel._points.size();

      for (iPoint = 0; iPoint < numberPoints; ++iPoint) {
      
	outputStream << "Point: " << iPoint << std::endl;
	outputStream << krigingModel._points[iPoint] << std::endl;
      
	std::vector<Value>::size_type iValue;
	const std::vector<Value>::size_type numberValues = 
	  krigingModel._values[iPoint].size();
      
	for (iValue = 0; iValue < numberValues; ++iValue) {

	  outputStream << "Value: " << std::endl;
	  outputStream << krigingModel._values[iPoint][iValue]
		       << std::endl;

	}
	
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

      //
      // firewalls
      //

      assert(krigingModel._AZ.size() == krigingModel._BZ.size());
      assert(krigingModel._AZ.size() == krigingModel._sigmaSqr.size());

      std::vector<Vector>::size_type iValue;
      const std::vector<Vector>::size_type numberValues = 
	krigingModel._AZ.size();

      outputStream << "AZ: " << std::endl;

      for (iValue = 0; iValue < numberValues; ++iValue) {
      
	outputStream << "Value: " << iValue << std::endl;
	outputStream << krigingModel._AZ[iValue] << std::endl;
      
      }

      outputStream << "BZ: " << std::endl;
    
      for (iValue = 0; iValue < numberValues; ++iValue) {
      
	outputStream << "Value: " << iValue << std::endl;
	outputStream << krigingModel._BZ[iValue] << std::endl;
      
      }
      
      //
      // sigmaSqr
      //
    
      outputStream << "sigmaSqr: " << std::endl;

      for (iValue = 0; iValue < numberValues; ++iValue) {
      
	outputStream << "Value: " << iValue << std::endl;
	outputStream << krigingModel._sigmaSqr[iValue] << std::endl;
      
      }
    

      //
      //
      //

      return outputStream;

    }

    //
    // polymorphic output
    //

    void 
    MultivariateKrigingModel::print(std::ostream & outputStream) const
    {
      
      outputStream << *this;
      return;

    }

}


