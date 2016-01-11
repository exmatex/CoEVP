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
// File:        KrigingDataBase.cc
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: Interpolation database using kriging interpolation.
//

#ifdef HAVE_PKG_libprof 
#ifdef HAVE_PROFILE_H
#include <Profile.h>
#endif // HAVE_PROFILE_H
#endif // HAVE_PKG_libprof

#include <fstream>

#include "KrigingDataBase.h"
#include <kriging/SecondMoment.h>
#include <base/ResponsePoint.h>
#include "ApproxNearestNeighbors.h"
#include <mtreedb/MTree.h>
#include <toolbox/database/HDFDatabase.h>
#include <murmur3/MurmurHash3.h>


#define STRING_DIGITS 16
#define MURMUR_SEED 42

#ifndef DEBUG
#  define DEBUG 0
#endif

using namespace krigalg;

namespace krigcpl {

uint128_t saved_model_key;

    namespace {


       uint128_t getKeyHash(const ResponsePoint& point)
       {
          int point_size = point.size();

          std::vector<double> data;
          data.resize(point_size);
          for (int i=0; i<point_size; ++i) {
             data[i] = point[i];
          }

	  uint128_t hash;
          MurmurHash3_x64_128(&data[0], point_size*sizeof(double)/sizeof(char), MURMUR_SEED, &hash);

          return hash;
       }

       void unpackKey(std::string& key, std::vector<double>& data)
       {
          char buffer[80];
          const char * ptr = key.c_str();
          unsigned numChar_total = 0;
          unsigned numChar;
          unsigned data_length = 0;

          data.clear();

          while (sscanf(ptr+numChar_total, "%[^:]%n", buffer, &numChar) != EOF) {
             numChar++;
             numChar_total += numChar;
             data_length++;
             data.resize(data_length);
             //             data[data_length-1] = atof(buffer);
             data[data_length-1] = strtod(buffer,NULL);
          }
       }

      //
      // local functions
      //

      //
      // given a point find closest CoKrigingModel 
      //

      std::pair<int, InterpolationModelPtr>
      findClosestCoKrigingModel(const ResponsePoint        & point,
                                ApproxNearestNeighbors     & ann,
                                krigalg::InterpolationModelFactoryPointer modelFactory,
                                ModelDatabase * modelDB,
				double                       maxQueryPointModelDistance)
      {

#ifdef HAVE_PKG_libprof

	ProfileBegin("findClosest");

#endif // HAVE_PKG_libprof

	//
	// query tree for the closest model
	//

        int k_neighbors = 1;

        std::vector<int> ids(k_neighbors);
        std::vector<uint128_t> keys(k_neighbors);
        std::vector<double> dists(k_neighbors);

        std::vector<double> x;
        x.resize(point.size());
        for (int i=0; i<point.size(); ++i) {
           x[i] = point[i];
        }
        ann.knn(x, k_neighbors, ids, keys, dists);

        bool found_neighbors = (ids.size() == k_neighbors);

	InterpolationModelPtr closestKrigingModel;
	int closestKrigingModelId;

        //
        // If a neighbor is found, compare the distance between the query point and the model
        // with the value of maxQueryPointModelDistance
        //

	if ( found_neighbors &&
             dists[0] <= maxQueryPointModelDistance ) {

           //
           // get handle to the located object
           //

           uint128_t model_key = keys[0];

           closestKrigingModel = modelDB->extract(model_key);

           closestKrigingModelId = ids[0];
        }
        else {
           closestKrigingModelId = -1;
        }

#ifdef HAVE_PKG_libprof

        ProfileEnd("findClosest");

#endif // HAVE_PKG_libprof

	return std::make_pair(closestKrigingModelId,
			      closestKrigingModel);
	
      }

      //
      // compute kriging error at a point; the issue here is that if the
      // kriging model contains a single point then the error will
      // naturally be computed as zero; if this is the case we will try to
      // estimate the error wrt a constant function
      //

      double
      compKrigingError(const InterpolationModel & krigingModel,
		       const Point              & queryPoint,
		       int                        valueId,
		       double                     meanErrorFactor)
      {

	const int numberPoints = krigingModel.getNumberPoints();

	//
	// firewalls
	//

	assert(numberPoints >= 1);

	//
	// compute min number of points to attempt meaningful
	// interpolation
	//

	const int minNumberPoints = krigingModel.hasGradient() ? 1 : 
	  2*(krigingModel.getPointDimension() + 1) - 1;

	//
	// compute the error if the kriging model contains a single point;
	// otherwise, simply return the kriging prediction
	//

	if (numberPoints <= minNumberPoints ) {

	  //
	  // get the kriging estimate at query point
	  //

	  const Value queryValue = krigingModel.interpolate(valueId,
							    queryPoint);
      
	  //
	  // get all points in the model; there should really only be
	  // ONE if we have gradient information
	  //
      
	  const std::vector<Point> & points = krigingModel.getPoints();
	  assert( krigingModel.hasGradient() ? (points.size() == 1) : true );

	  //
	  // get the kriging estimate at the origin of the kriging model
	  //

	  const Value originValue = 
	    krigingModel.interpolate(valueId,
				     points.front());

	  //
	  // compute the error as the difference between the value at the
	  // queryPoint and origin point
	  //

	  return (queryValue[0] - originValue[0])*
	    (queryValue[0] - originValue[0]);


	} else 
	  return meanErrorFactor*meanErrorFactor*
	    krigingModel.getMeanSquaredError(valueId, queryPoint)[0];

	//
	// can never be reached
	//
    
	return 0.0;
    
      }

      //
      // use a kriging model to compute the infinity norm of the
      // mean squared error
      //

      double
      checkErrorInf(InterpolationModelPtr   krigingModel,
		    const ResponsePoint   & queryPoint,
		    int                     valueDimension,
		    double                 _meanErrorFactor)
      {

	double maxError = 0.0;

	//
	// iterate over values
	//

	for (int iValue = 0; iValue < valueDimension; ++iValue) {
	  
	  //
	  // compute the error estimate
	  //
		  
	  const double errorEstimate = 
	    compKrigingError(*krigingModel, 
			     queryPoint,
			     iValue,
			     _meanErrorFactor);
	  
	  //
	  // save error
	  //

	  maxError = std::max(maxError,
			      errorEstimate);


	}

	return maxError;

      }

      double
      checkErrorL2(InterpolationModelPtr   krigingModel,
		   const ResponsePoint   & queryPoint,
		   int                     valueDimension,
		   double                 _meanErrorFactor)
      {

	double maxError = 0.0;

	//
	// iterate over values
	//

	for (int iValue = 0; iValue < valueDimension; ++iValue) {
	  
	  //
	  // compute the error estimate
	  //
		  
	  const double errorEstimate = 
	    compKrigingError(*krigingModel, 
			     queryPoint,
			     iValue,
			     _meanErrorFactor);
	  
	  //
	  // save error
	  //

	  maxError += errorEstimate;

	}

	return maxError/valueDimension;

      }

      inline double
      checkError(InterpolationModelPtr   krigingModel,
		 const ResponsePoint   & queryPoint,
		 int                     valueDimension,
		 double                 _meanErrorFactor)
      {

	return checkErrorInf(krigingModel,
			     queryPoint,
			     valueDimension,
			     _meanErrorFactor);

      }


      //
      // Given a point find the "best" kriging model. Here, best means
      // that its use results in the smallest possible interpolation
      // error. This function returns a handle to the best model that
      // satisfies given tolerance or id_undefined if no model can be found.
      //

      std::pair<int, InterpolationModelPtr>
      findBestCoKrigingModel(bool &                       canInterpolateFlag,
			     const ResponsePoint &        point,
                             ApproxNearestNeighbors     & ann,
                             ModelDatabase * modelDB,
                             const InterpolationModelFactoryPointer& _modelFactory,
			     double                       tolerance,
			     double                       meanErrorFactor,
			     double                       maxQueryPointModelDistance,
			     int                          maxNumberSearchModels,
			     int                          valueDimension)
      {
	
	canInterpolateFlag = false;

	//
	// query the tree for the maxNumberSearchModels closest models
	//

        int k_neighbors = maxNumberSearchModels;

        std::vector<int> ids(k_neighbors);
        std::vector<uint128_t> keys(k_neighbors);
        std::vector<double> dists(k_neighbors);

        std::vector<double> x;
        x.resize(point.size());
        for (int i=0; i<point.size(); ++i) {
           x[i] = point[i];
        }
        ann.knn(x, k_neighbors, ids, keys, dists);

        bool found_neighbors = (ids.size() > 0);

	InterpolationModelPtr bestKrigingModel;
	int bestKrigingModelId;

        std::pair<int, InterpolationModelPtr> bestCoKrigingModel;

        if ( found_neighbors ) {

           //
           // iterate through the search results
           //

           double minError = std::numeric_limits<double>::max();

           for (int iter=0; iter<ids.size(); ++iter) {

              uint128_t model_key = keys[iter];

              InterpolationModelPtr krigingModel = modelDB->extract(model_key);
              //
              // skip invalid models
              //

              if (krigingModel->isValid() == false)
                 continue;

              //
              // compute error
              //

              const double errorEstimate = checkError(krigingModel,
                                                      point,
                                                      valueDimension,
                                                      meanErrorFactor);
              
              //
              // check if a better model encountered
              //
              
              if (errorEstimate < minError) {
                 
                 bestCoKrigingModel = std::make_pair(ids[iter], krigingModel);
                 
                 minError = errorEstimate;
                 
              }
              
           }
	 
           //
           // check if the best model satisfies tolerance requirement
           //
	
           if (minError <= tolerance*tolerance)
              canInterpolateFlag = true;

        }
        else {
	  bestCoKrigingModel = make_pair(-1, InterpolationModelPtr());
        }

        return bestCoKrigingModel;

      }

      std::pair<int, InterpolationModelPtr>
      findBestCoKrigingModel(bool &                       canInterpolateFlag,
			     const ResponsePoint &        point,
                             ApproxNearestNeighbors     & ann,
                             ModelDatabase * modelDB,
                             const InterpolationModelFactoryPointer& _modelFactory,
			     double                       tolerance,
			     double                       meanErrorFactor,
			     double                       maxQueryPointModelDistance,
			     int                          maxNumberSearchModels,
			     int                          maxKrigingModelSize,
			     int                          valueDimension)
      {

#ifdef HAVE_PKG_libprof

	ProfileBegin("findBest");

#endif // HAVE_PKG_libprof

	canInterpolateFlag = false;

#ifdef HAVE_PKG_libprof

	  ProfileEnd("findBest");
	
#endif // HAVE_PKG_libprof

	//
	// query the tree for the maxNumberSearchModels closest models
	//

        int k_neighbors = maxNumberSearchModels;

        std::vector<int> ids(k_neighbors);
        std::vector<uint128_t> keys(k_neighbors);
        std::vector<double> dists(k_neighbors);

        std::vector<double> x;
        x.resize(point.size());
        for (int i=0; i<point.size(); ++i) {
           x[i] = point[i];
        }
        ann.knn(x, k_neighbors, ids, keys, dists);

        bool found_neighbors = (ids.size() > 0);

        if ( found_neighbors ) {

           //
           // iterate through the search results
           //

           for (int iter=0; iter<ids.size(); ++iter) {
              
              uint128_t model_key = keys[iter];

              InterpolationModelPtr krigingModel = modelDB->extract(model_key, (InterpolationModelFactoryPointer *)&_modelFactory);

              //
              // skip if invalid model
              //
              
              if (krigingModel->isValid() == false)
                 continue;

              //
              // skip if too far away; as results come back closest-first, we know
              // all of the other results are too far away as well, so can break
              // could probably 
              //
              
              if (dists[iter] > maxQueryPointModelDistance)
                 break;

              //
              // compute error
              //

              const double errorEstimate = checkError(krigingModel,
                                                      point,
                                                      valueDimension,
                                                      meanErrorFactor);

              if (errorEstimate <= tolerance*tolerance) {

#ifdef HAVE_PKG_libprof

                 ProfileEnd("findBest");
	
#endif // HAVE_PKG_libprof

                 canInterpolateFlag = true;

                 return std::make_pair(ids[iter], krigingModel);

              }
           }

           //
           // a model suitable for interpolation has not been
           // found-return the closest model
           //

           uint128_t model_key = keys[0];

           InterpolationModelPtr krigingModel = modelDB->extract(model_key, (InterpolationModelFactoryPointer *)&_modelFactory);

#ifdef HAVE_PKG_libprof

           ProfileEnd("findBest");
	
#endif // HAVE_PKG_libprof

           return std::make_pair(ids[0], krigingModel);

        }
        else {
           return std::make_pair(-1, InterpolationModelPtr());
        }
      }


      //
      // check error and interpolate using given kriging model
      //

      bool
      checkErrorAndInterpolate(double               * value,
			       InterpolationModelPtr  krigingModel,
			       const ResponsePoint  & queryPoint,
			       int                    valueDimension,
			       double                _tolerance,
			       double                _meanErrorFactor,
                               double               & errorEstimate)
      {
	const double toleranceSqr = _tolerance*_tolerance;
	
        errorEstimate = 0.;

	for (int iValue = 0; iValue < valueDimension; ++iValue) {
	  
	  //
	  // compute the error estimate
	  //
	
	  
#ifdef HAVE_PKG_libprof
	  
	  ProfileBegin("errEst");
	  
#endif // HAVE_PKG_libprof
	  
	  const double iErrorEstimate = 
	    compKrigingError(*krigingModel, 
			     queryPoint,
			     iValue,
			     _meanErrorFactor);
	  
          errorEstimate = std::max(errorEstimate, sqrt(iErrorEstimate));


#ifdef HAVE_PKG_libprof
	  
	  ProfileEnd("errEst");
	  
#endif // HAVE_PKG_libprof
	  
	  //
	  // check the errorEstimate against the tolerance; if the 
	  // estimate is greater than tolerance simpoy return failure
	  //
	  
          //	  if ( iErrorEstimate > toleranceSqr) {
	  if ( fabs(iErrorEstimate) > toleranceSqr) {

	    return false;
	    
	  }

#ifdef HAVE_PKG_libprof

	  ProfileBegin("interpolate");
	  
#endif // HAVE_PKG_libprof
	  
	  //
	  // compute the value
	  //
	  
	  const Value valueEstimate = 
	    krigingModel->interpolate(iValue,
				      queryPoint);
	  
	  //
	  // put the value of the function into value (valueEstimate
	  // contains the value of the function follwed by the gradient;
	  // here we are interested only in the value of the function so
	  // the gradient is simply discarded).
	  //

	  value[iValue] = valueEstimate[0];
	  
	  // 	std::cout << "value id: " << iValue 
	  // 		  << " error estimate: " << sqrt(iErrorEstimate) 
	  // 		  << std::endl;
          //          std::cout << "value id: " << iValue 
          //                    << " value estimate: " << valueEstimate[0]
          //                    << std::endl;

#ifdef HAVE_PKG_libprof
	  
	  ProfileEnd("interpolate");
	  
#endif // HAVE_PKG_libprof
	  
	}
	
	return true;

      }
			       
      bool
      checkErrorAndInterpolate(double               * value,
			       double               * gradient,
			       InterpolationModelPtr  krigingModel,
			       const ResponsePoint  & queryPoint,
			       int                    pointDimension,
			       int                    valueDimension,
			       double                _tolerance,
			       double                _meanErrorFactor,
                               double               & errorEstimate )
      {

	//
	// firewall
	//

	assert(krigingModel->hasGradient() == true);

	//
	// estimate error for all values from the kriging models
	//

	const double toleranceSqr = _tolerance*_tolerance;

        errorEstimate = 0.;

	for (int iValue = 0; iValue < valueDimension; ++iValue) {

#ifdef HAVE_PKG_libprof

	  ProfileBegin("errEst");
	
#endif // HAVE_PKG_libprof
	
	  //
	  // compute the error estimate
	  //
	
	  const double iErrorEstimate = 
	    compKrigingError(*krigingModel, 
			     queryPoint,
			     iValue,
			     _meanErrorFactor);

          errorEstimate = std::max(errorEstimate, sqrt(iErrorEstimate));

#ifdef HAVE_PKG_libprof

	  ProfileEnd("errEst");

#endif // HAVE_PKG_libprof

	  //
	  // check the errorEstimate against the tolerance; if the 
	  // estimate is greater than tolerance simpoy return failure
	  //

          //	  if ( iErrorEstimate > toleranceSqr) {
	  if ( fabs(iErrorEstimate) > toleranceSqr) {

	    return false;

	  }

#ifdef HAVE_PKG_libprof

	  ProfileBegin("interpolate");

#endif // HAVE_PKG_libprof

	  //
	  // compute the value
	  //

	  const Value valueEstimate = 
	    krigingModel->interpolate(iValue,
				      queryPoint);

	  //
	  // put the value of the function into value (valueEstimate
	  // contains the value of the function follwed by the gradient;
	  //

	  value[iValue] = valueEstimate[0];

	  //
	  // store gradient data
	  //
	  
	  for (int i = 0; i < pointDimension; ++i)
	    gradient[i*valueDimension + iValue] = valueEstimate[1 + i];

#ifdef HAVE_PKG_libprof
	  
	  ProfileEnd("interpolate");

#endif // HAVE_PKG_libprof

	}

	return true;

      }

      //
      // Use a given kriging model to interpolate. The assumption is
      // that the kriging model is the "best" possible choice so there
      // is no reason to check that the error is withing given tolerance 
      // as above (checkErrorAndInterpolate()).
      //

      void
      interpolate(double                * value,
		  InterpolationModelPtr   krigingModel,
		  const ResponsePoint   & queryPoint,
		  int                     valueDimension)
      {

#ifdef HAVE_PKG_libprof

	ProfileBegin("best_interp");

#endif // HAVE_PKG_libprof

	for (int iValue = 0; iValue < valueDimension; ++iValue) {

	  //
	  // compute the value
	  //
	  
	  const Value valueEstimate = 
	    krigingModel->interpolate(iValue,
				      queryPoint);
	  
	  //
	  // put the value of the function into value (valueEstimate
	  // contains the value of the function follwed by the gradient;
	  // here we are interested only in the value of the function so
	  // the gradient is simply discarded).
	  //

	  value[iValue] = valueEstimate[0];
	  
	}

#ifdef HAVE_PKG_libprof

	ProfileEnd("best_interp");

#endif // HAVE_PKG_libprof


	return;

      }

      void
      interpolate(double                                * value,
		  double                                * gradient,
		  InterpolationModelPtr                   krigingModel,
		  const ResponsePoint                   & queryPoint,
		  int                                     pointDimension,
		  int                                     valueDimension)
      {

#ifdef HAVE_PKG_libprof

	ProfileBegin("best_interp");

#endif // HAVE_PKG_libprof

	//
	// firewalls
	//

	assert(krigingModel->hasGradient() == true);

	for (int iValue = 0; iValue < valueDimension; ++iValue) {

	  //
	  // compute the value
	  //

	  const Value valueEstimate = 
	    krigingModel->interpolate(iValue,
				      queryPoint);

	  //
	  // put the value of the function into value (valueEstimate
	  // contains the value of the function follwed by the gradient;
	  //

	  value[iValue] = valueEstimate[0];

	  //
	  // store gradient data
	  //
	  
	  for (int i = 0; i < pointDimension; ++i)
	    gradient[i*valueDimension + iValue] = valueEstimate[1 + i];

	}

#ifdef HAVE_PKG_libprof

	ProfileEnd("best_interp");

#endif // HAVE_PKG_libprof


	return;

      }

      //
      // convert value data to make it ready for kriging model
      // insertion
      //

      std::vector<Value> 
      copyValueData(const double * value,
		    const double * gradient,
		    int            pointDimension,
		    int            valueDimension)
      {

	//
	// instatiate return object
	//

	std::vector<Value> pointValues;

	//
	// iterate over all values
	//

	for (int iValue = 0; iValue < valueDimension; ++iValue) {

	  //
	  // copy value data
	  //

	  Value pointValue(pointDimension + 1);
	  
	  pointValue[0] = value[iValue];

	  //
	  // copy gradient data
	  //

	  for (int i = 0; i < pointDimension; ++i)
	    pointValue[1 + i] = gradient[i*valueDimension + iValue];

	  //
	  // add pointValue
	  //

	  pointValues.push_back(pointValue);

	}

	return pointValues;

      }

      //
      // convert value data to make it ready for kriging model
      // insertion (no gradient version)
      //

      std::vector<Value> 
      copyValueData(const double * value,
		    int            pointDimension,
		    int            valueDimension)
      {

	//
	// instatiate return object
	//

	std::vector<Value> pointValues;

	//
	// iterate over all values
	//

	for (int iValue = 0; iValue < valueDimension; ++iValue) {

	  //
	  // copy value data
	  //

	  Value pointValue(1);
	  
	  pointValue[0] = value[iValue];

	  //
	  // add pointValue
	  //

	  pointValues.push_back(pointValue);

	}

	return pointValues;

      }

      //
      // instantiate, fill in and insert a new kriging model
      //

       void
       addNewModel(
                   ModelDatabase * modelDB,
                   ApproxNearestNeighbors&                  ann,
                   const InterpolationModelFactoryPointer & _modelFactory,
                   int &                                    objectId,
                   const double *                           pointData,
                   const double *                           valueData,
                   const double *                           gradientData,
                   int                                      pointDimension,
                   int                                      valueDimension)
      {
	
	//
	// start new kriging model
	//
	
	InterpolationModelPtr krigingModel = _modelFactory->build();

	//
	// copy value and gradient data
	//
	
	std::vector<Value> pointValue;


	if (krigingModel->hasGradient() == true) 
	  pointValue = copyValueData(valueData,
				     gradientData,
				     pointDimension,
				     valueDimension);
	else
	  pointValue = copyValueData(valueData,
				     pointDimension,
				     valueDimension);

	//
	// copy point data into a Point object
	//
	
	const ResponsePoint point(pointDimension,
				  pointData );
	
	//
	// insert point/value data into the model
	//
	
	krigingModel->addPoint(point,
			       pointValue);

        // Create a key string corresponding to the new point at
        // which the new interpolation model is centered

        // Insert the model key into the approximate nearest neighbor database

        std::vector<double> point_data;
        point_data.resize(point.size());
        for (int i=0; i<point.size(); ++i) {
           point_data[i] = point[i];
        }

        uint128_t model_key = getKeyHash(point);
        objectId = ann.insert(point_data, model_key);

        // Insert the interpolation model into the interpolation model database

        modelDB->insert(model_key,krigingModel, (ResponsePoint *)&point);

	return;
	
      }

      //
      // compute center of mass for a kriging model
      //

      Point
      getModelCenterMass(const InterpolationModel & krigingModel)
      {

	//
	// get all points in the model
	//

	const std::vector<Point> & points = krigingModel.getPoints();
    
	//
	// firewall
	//

	assert(points.empty() == false);

	//
	// instantiate center od mass
	//

	Point centerMass(points.front().size(),
			 0.0);

	//
	// iterate over all points
	//
    
	std::vector<Point>::const_iterator pointsIter;
	std::vector<Point>::const_iterator pointsEnd = points.end();
    
	for (pointsIter  = points.begin();
	     pointsIter != pointsEnd;
	     ++pointsIter) {
      
	  //
	  // get Point handle
	  //
      
	  const Point & point = *pointsIter;
      
	  //
	  // accumulate 
	  //
      
	  centerMass += static_cast<Vector>(point);

	}
    
	//
	// scale
	//
    
	mtl::scale(centerMass, 
		   1.0/points.size());
	
	return centerMass;

      }

      //
      // print kriging point statistics
      //

#if 0
      inline void
      outputKrigingModelStats(std::ostream & outputStream,
			      const DB  &    keyDB,
			      int            maxKrigingModelSize)
      {
	
	return;

      }
#endif

      //
      // read the contents of the data store and insert into the tree
      //

      std::pair<int, int>
      initializeModelDBFromFile(
                                ApproxNearestNeighbors&                  ann,
				const InterpolationModelFactoryPointer & _modelFactory,
				const std::string &                      directoryName,
				const std::string &                      prefix)
      {
#ifdef HAVE_PKG_hdf5
	//
	// open and mount data store summary
	//

	toolbox::HDFDatabase summaryDatabase("data_store_summary_database");
	
	summaryDatabase.mount(directoryName + "/" + prefix +
			      "__data_store_summary",
			      "R");

	//
	// get number of objects
	//

	const int numberObjects = summaryDatabase.getInteger("d_num_objects");
	const int objectFileCapacity = 
	  summaryDatabase.getInteger("d_object_file_capacity");

	//
	// aggregate of all point/value pairs in the kriging model
	//

	int totalNumberPoints = 0;

	//
	// iterate over all files and read kriging models
	//
	
	const int numberFiles = numberObjects/objectFileCapacity + 1;

	for (int iFileIndex = 0; iFileIndex < numberFiles; ++iFileIndex) {

	  //
	  // construct file name
	  //

	  std::stringstream fileStringStream;

	  fileStringStream << directoryName << "/" << prefix 
			   << "__data_objects."
			   << std::setfill('0') << std::setw(8) 
			   << iFileIndex;
	  
	  
	  std::string iFileName;
	  
	  fileStringStream >> iFileName;
	  
	  //
	  // open and mount model database
	  //
	  
	  toolbox::HDFDatabase modelDatabase("model_database");
	  
	  modelDatabase.mount(iFileName,
			      "R");
	  
	  //
	  // determine the number of models in the file
	  //
	  
	  const int numberModelsInFile = 
	    (iFileIndex == numberFiles - 1) ? 
	    numberObjects - iFileIndex*objectFileCapacity : 
	    objectFileCapacity;

	  //
	  // iterate all objects in the database
	  //

	  for (int iModel = 0; iModel < numberModelsInFile; ++iModel) {

	    //
	    // determine model id
	    //

	    const int modelId = iFileIndex*objectFileCapacity + iModel;

	    //
	    // construct object string 
	    //

	    std::stringstream modelStringStream;

	    modelStringStream << "data_object__" 
			      << std::setfill('0') << std::setw(9)
			      << modelId;

	    std::string modelString;
	    
	    modelStringStream >> modelString;
	    
	    //
	    // mount modelString database
	    //

	    toolbox::DatabasePtr objectDatabase = 
	      modelDatabase.getDatabase(modelString);
	      
	    //
	    // instantiate a new model
	    //

	    krigalg::InterpolationModelPtr krigingModelPtr = 
	      _modelFactory->build();

	    //
	    // fill the contents of the new model
	    //

	    krigingModelPtr->getFromDatabase(*objectDatabase);

	    //
	    // update totalNumberPoints
	    //

	    totalNumberPoints += krigingModelPtr->getNumberPoints();

	    //
	    // get the center of the kriging model
	    //

	    Point krigingModelCenter = getModelCenterMass(*krigingModelPtr);

	    //
	    // copy krigingModelCenter into ResponsePoint
	    //

	    const ResponsePoint point(krigingModelCenter.size(),
				      &(krigingModelCenter[0]));

            // Create a key corresponding to interpolation model center
            // and insert it into the approximate nearest neighbor database

            std::vector<double> point_data;
            point_data.resize(point.size());
            for (int i=0; i<point.size(); ++i) {
               point_data[i] = point[i];
            }

            uint128_t model_key = getKeyHash(point);
            ann.insert(point_data, model_key);

            // Insert the interpolation model into the interpolation model database

            this->_modelDB->insert(model_key, krigingModelPtr, (ResponsePoint *)&point);

	  }

	  //
	  // unmount modelDatabase
	  //

	  modelDatabase.unmount();

	}

	//
	// cleanup
	//

	summaryDatabase.unmount();

	return std::make_pair(numberObjects, totalNumberPoints);

#else
	TBOX_ERROR("KrigingInterpolationDatabase error..."
		   << "\n   Cannot read model database from file..."
		   << " code not compiled with HDF" << endl);

	return std::make_pair(-1,-1);
#endif
      }

      //
      // output kriging models position data
      //

#if 0
      void
      outputKrigingModelPositionData(const std::string & fileName,
				     DB                & keyDB)
      {

         if (typeid(keyDB) == typeid(MTree)) {

            //
            // make sure that the tree has been properly intialized; at
            // the moment the only way to do this is to check the number
            // of levels
            //

            if (((MTree&)keyDB).getNumberLevels() == 0)
               return;

            //
            // open output stream
            //

            std::ofstream outputStream(fileName.c_str());
	
            //
            // we rely here on MTree stats
            //

            ((MTree&)keyDB).calculateLevelStatistics();

            //
            // get stats for level 0 (leaf nodes)
            //

            const MTreeLevelStatistic * mtreeStats =
               ((MTree&)keyDB).getLevelStatistics(0);

            //
            // get number of leaf nodes
            //

            const int numberLeafNodes = mtreeStats->getNumberNodesOnLevel();
	
            //
            // iterate over leaf nodes
            //

            for (int iNode = 0; iNode < numberLeafNodes; ++iNode) {
	  
               //
               // get handle to node stat
               //

               const MTreeNodeStat & leafNodeStat = mtreeStats->getNodeStat(iNode);
	 
               //
               // get handle to data objects
               //

               const std::vector<int> & krigingModelObjectIds = 
                  leafNodeStat.getDataObjectIds();

               //
               // iterate over data object Ids
               //

               std::vector<int>::const_iterator krigingModelObjectIdsEnd = 
                  krigingModelObjectIds.end();
               std::vector<int>::const_iterator krigingModelObjectIdsIter;

               for (krigingModelObjectIdsIter  = krigingModelObjectIds.begin();
                    krigingModelObjectIdsIter != krigingModelObjectIdsEnd;
                    ++krigingModelObjectIdsIter) {

                  //
                  // get pointer to the MTreeObject
                  //

                  const DBObjectPtr dbObjectPtr = 
                     ((MTree&)keyDB).getObject(*krigingModelObjectIdsIter);

                  //
                  // get handle to the underlying DBKrigingModelObject
                  //

                  const DBKrigingModelObject & dbObject = 
                     dynamic_cast<const DBKrigingModelObject &>(*dbObjectPtr);

                  //
                  // get pointer to the kriging model
                  //

                  const InterpolationModelPtr krigingModel = dbObject.getModel();

                  //
                  // compute center of the kriging model
                  //

                  const Point modelCenter = getModelCenterMass(*krigingModel);

                  //
                  // output model center
                  //
	    
                  outputStream << modelCenter;
                  
               }
            }
         }
         else {

            // Need to figure out how to generalize this to work
            // with databases other than MTree

         }

         return;
	
      }
#endif
    }

    // Object class member definitions start here

    //
    // construction/destruction
    //

    KrigingDataBase::KrigingDataBase(int pointDimension,
                                     int valueDimension,
                                     const InterpolationModelFactoryPointer  & modelFactory,
                                     ApproxNearestNeighbors& ann,
                                     int    maxKrigingModelSize,
                                     int    maxNumberSearchModels,
                                     bool   useHint,
                                     double meanErrorFactor,
                                     double tolerance,
                                     double maxQueryPointModelDistance,
                                     int    agingThreshold, 
                                     ModelDatabase * modelDB)
      : InterpolationDataBase(pointDimension,
			      valueDimension),
	_modelFactory(modelFactory),
	_maxKrigingModelSize(maxKrigingModelSize),
	_maxNumberSearchModels(maxNumberSearchModels),
	_useHint(useHint),
	_meanErrorFactor(meanErrorFactor),
	_tolerance(tolerance),
	_maxQueryPointModelDistance(maxQueryPointModelDistance),
        _ann(ann),
	_numberKrigingModels(0),
	_numberPointValuePairs(0),
	_agingThreshold(agingThreshold),
	_modelDB(modelDB)
    {
      return;
      
    }

    KrigingDataBase::KrigingDataBase(int pointDimension,
                                     int valueDimension,
                                     const InterpolationModelFactoryPointer  & modelFactory,
                                     ApproxNearestNeighbors& ann,
                                     int    maxKrigingModelSize,
                                     int    maxNumberSearchModels,
                                     bool   useHint,
                                     double meanErrorFactor,
                                     double tolerance,
                                     double maxQueryPointModelDistance,
                                     int    agingThreshold,
                                     const std::string & directoryName,
                                     const std::string & fileName,
                                     ModelDatabase * modelDB)
       : InterpolationDataBase(pointDimension,
                               valueDimension,
                               fileName),
         _modelFactory(modelFactory),
         _maxKrigingModelSize(maxKrigingModelSize),
         _maxNumberSearchModels(maxNumberSearchModels),
         _useHint(useHint),
         _meanErrorFactor(meanErrorFactor),
         _tolerance(tolerance),
         _maxQueryPointModelDistance(maxQueryPointModelDistance),
         _ann(ann),
         _numberKrigingModels(0),
         _numberPointValuePairs(0),
         _agingThreshold(agingThreshold),
         _modelDB(modelDB)
    {
      return;
       const std::pair<int, int> kriginigModelsStats =
          initializeModelDBFromFile(
                                    _ann,
                                    _modelFactory,
                                    directoryName,
                                    fileName);

       _numberKrigingModels   += kriginigModelsStats.first;
       _numberPointValuePairs += kriginigModelsStats.second;



       return;
      
    }

    KrigingDataBase::~KrigingDataBase()
    {

      return;

    }
    
    //
    // Compute interpolated value at a point.
    //

    // Version 1 (no derivatives)

    bool
    KrigingDataBase::interpolate(double            * value,
                                 int               & hint,
                                 const double      * point,
                                 std::vector<bool> & flags,
                                 double            & error_estimate)
    {

      //
      // make sure there is enough space in flags 
      //

      assert(flags.size() >= NUMBER_FLAGS);

      //
      // initialize flags container
      //

      std::fill(flags.begin(),
		flags.end(),
		false);

      //
      // shortcuts to frequently accesses data
      //

      const int pointDimension = getPointDimension();
      const int valueDimension = getValueDimension();

      //
      // instatiate point object from point data
      //

      const ResponsePoint queryPoint(pointDimension,
				     point);

      //
      // use hint to get the most recently used model
      //

      if (_useHint && hint != -1) {
	
         uint128_t model_key = _ann.getKey(hint);

         if (model_key == uint128_t_undefined) {

            flags[LOST_HINT_FLAG] = true;

         } else {
            InterpolationModelPtr hintKrigingModel = this->_modelDB->extract(model_key, (InterpolationModelFactoryPointer *)&this->_modelFactory);

            //
            // check if can interpolate; need a valid model for this
            //

            if (hintKrigingModel->isValid() == true) {

               //
               // check the distance between hintKrigingModel and point
               //
               const Point modelCenter = getModelCenterMass(*hintKrigingModel);
               const Vector pointRelativePosition = queryPoint - modelCenter;
               const double distanceSqr = krigalg::dot(pointRelativePosition,
                                                       pointRelativePosition);

               if (distanceSqr > 
                   _maxQueryPointModelDistance*_maxQueryPointModelDistance) {
                  flags[LOST_HINT_FLAG] = true;
               } else {

                  const bool hintModelSuccess = 
                     checkErrorAndInterpolate(value,
                                              hintKrigingModel,
                                              queryPoint,
                                              valueDimension,
                                              _tolerance,
                                              _meanErrorFactor,
                                              error_estimate );

                  if (hintModelSuccess == true) {
                     flags[USED_HINT_FLAG] = true;
                     return true;
                  }
	      
               }
               
            }
            
         }
	
      }

      //
      // A kriging model based on hint did not produce a valid interpolation.
      // Find closest kriging model.
      //

      if (_maxNumberSearchModels == 1) {

#ifdef HAVE_PKG_libprof

	ProfileBegin("interpClosest");

#endif // HAVE_PKG_libprof

	const std::pair<int, InterpolationModelPtr> 
	  closestKrigingModelData = findClosestCoKrigingModel(queryPoint,
                                                              _ann,
                                                              _modelFactory,
                                                              _modelDB,
							      _maxQueryPointModelDistance);
      
	InterpolationModelPtr closestKrigingModel = 
	  closestKrigingModelData.second;
	hint = closestKrigingModelData.first;

	//
	// if no kriging model is available return
	//

	if (hint == id_undefined || closestKrigingModel->isValid() == false) {

#ifdef HAVE_PKG_libprof

	  ProfileEnd("interpClosest");

#endif // HAVE_PKG_libprof

	  return false;

	}

	//
	// estimate error for all values from the kriging models
	//
      
	const bool interpolationSuccess =  
	  checkErrorAndInterpolate(value,
				   closestKrigingModel,
				   queryPoint,
				   valueDimension,
				   _tolerance,
				   _meanErrorFactor,
                                   error_estimate);

#ifdef HAVE_PKG_libprof

	ProfileEnd("interpClosest");

#endif // HAVE_PKG_libprof

	return interpolationSuccess;

      } else {

#ifdef HAVE_PKG_libprof

	ProfileBegin("interpBest");

#endif // HAVE_PKG_libprof

	//
	// search more than one kriging model
	//

	bool canInterpolateFlag;

	const std::pair<int, InterpolationModelPtr> 
	  bestKrigingModelData = findBestCoKrigingModel(canInterpolateFlag,
							queryPoint,
                                                        _ann,
                                                        _modelDB,
                                                        _modelFactory,
							_tolerance,
							_meanErrorFactor,
							_maxQueryPointModelDistance,
							_maxNumberSearchModels,
							_maxKrigingModelSize,
							valueDimension);
	
	InterpolationModelPtr bestKrigingModel = 
	  bestKrigingModelData.second;
	hint = bestKrigingModelData.first;

	//
	// if no kriging model is available return
	//

	if (canInterpolateFlag == false) {

#ifdef HAVE_PKG_libprof

	  ProfileEnd("interpBest");

#endif // HAVE_PKG_libprof

	  return false;

	}

	//
	// interpolate using the best kriging model available
	//

 	return checkErrorAndInterpolate(value,
 					bestKrigingModel,
 					queryPoint,
 					valueDimension,
 					_tolerance,
 					_meanErrorFactor,
                                        error_estimate);

#ifdef HAVE_PKG_libprof

	ProfileEnd("interpBest");

#endif // HAVE_PKG_libprof

	return true;
	
      }

      //
      // this should never be reached
      //

      assert(false);

      return false;

    }

    // Version 2

    bool 
    KrigingDataBase::interpolate(double            * value,
                                 double            * gradient,
                                 int               & hint,
                                 const double      * point,
                                 std::vector<bool> & flags,
                                 double            & error_estimate)
    {

      //
      // make sure there is enough space in flags 
      //

      assert(flags.size() >= NUMBER_FLAGS);

      //
      // initialize flags container
      //

      std::fill(flags.begin(),
		flags.end(),
		false);
      //
      // shortcuts to frequently accesses data
      //

      const int pointDimension = getPointDimension();
      const int valueDimension = getValueDimension();

      //
      // instatiate point object from point data
      //

      const ResponsePoint queryPoint(pointDimension,
				     point);

      //
      // try hint (if valid)
      //

      if (_useHint && hint != -1) {
	
         uint128_t model_key = _ann.getKey(hint);

        if (model_key == uint128_t_undefined) {

	  flags[LOST_HINT_FLAG] = true;

	} else {
    
           InterpolationModelPtr hintKrigingModel = _modelDB->extract(model_key, (InterpolationModelFactoryPointer *)&this->_modelFactory);
	  
	  //
	  // check the distance between hintKrigingModel and point
	  //
	  const Point modelCenter = getModelCenterMass(*hintKrigingModel);
	  const Vector pointRelativePosition = queryPoint - modelCenter;
	  const double distanceSqr = krigalg::dot(pointRelativePosition,
	  					    pointRelativePosition);

	  if (distanceSqr > 
	  	_maxQueryPointModelDistance*_maxQueryPointModelDistance) {
	    flags[LOST_HINT_FLAG] = true;
	  } else {
	  
	    const bool hintModelSuccess = 
	      checkErrorAndInterpolate(value,
                                       gradient,
                                       hintKrigingModel,
                                       queryPoint,
                                       pointDimension,
                                       valueDimension,
                                       _tolerance,
                                       _meanErrorFactor,
                                       error_estimate);
	    
// 	    if (hintModelSuccess == true)
// 	      std::cout << "hint success" << std::endl;
	    
	    if (hintModelSuccess == true) {
	      flags[USED_HINT_FLAG] = true;
	      return true;
	    }

	  }

	}
	
      }
      
      //
      // find closest kriging model
      //

      if (_maxNumberSearchModels == 1) {

#ifdef HAVE_PKG_libprof

	ProfileBegin("interpClosest");

#endif // HAVE_PKG_libprof

	const std::pair<int, InterpolationModelPtr> 
	  closestKrigingModelData = findClosestCoKrigingModel(queryPoint,
                                                              _ann,
                                                              _modelFactory,
                                                              _modelDB,
							      _maxQueryPointModelDistance);
      
	InterpolationModelPtr closestKrigingModel = 
	  closestKrigingModelData.second;
	hint = closestKrigingModelData.first;

	//
	// if no kriging model is available return
	//

	if (hint == id_undefined) {

#ifdef HAVE_PKG_libprof

	ProfileEnd("interpClosest");

#endif // HAVE_PKG_libprof

	  return false;

	}
	  
	//
	// check error and interpolate
	//

	const bool interpolationSuccess = 
	  checkErrorAndInterpolate(value,
				   gradient,
				   closestKrigingModel,
				   queryPoint,
				   pointDimension,
				   valueDimension,
				   _tolerance,
				   _meanErrorFactor,
                                   error_estimate);

#ifdef HAVE_PKG_libprof

	ProfileEnd("interpClosest");

#endif // HAVE_PKG_libprof

	return interpolationSuccess;

      } else {

#ifdef HAVE_PKG_libprof

	ProfileBegin("interpBest");

#endif // HAVE_PKG_libprof

	//
	// search more than one kriging model
	//

	bool canInterpolateFlag;

	const std::pair<int, InterpolationModelPtr> 
	  bestKrigingModelData = findBestCoKrigingModel(canInterpolateFlag,
							queryPoint,
                                                        _ann,
                                                        _modelDB,
                                                        _modelFactory,
							_tolerance,
							_meanErrorFactor,
							_maxQueryPointModelDistance,
							_maxNumberSearchModels,
							_maxKrigingModelSize,
							valueDimension);
	
	InterpolationModelPtr bestKrigingModel = 
	  bestKrigingModelData.second;
	hint = bestKrigingModelData.first;

	//
	// if no kriging model is available return
	//

	if (canInterpolateFlag == false) {

#ifdef HAVE_PKG_libprof

	  ProfileEnd("interpBest");

#endif // HAVE_PKG_libprof

	  return false;

	}

	//
	// interpolate using the best kriging model available
	//

        krigcpl::interpolate(value,
                             gradient,
                             bestKrigingModel,
                             queryPoint,
                             pointDimension,
                             valueDimension);

        // 	return checkErrorAndInterpolate(value,
        // 					gradient,
        // 					bestKrigingModel,
        // 					queryPoint,
        // 					pointDimension,
        // 					valueDimension,
        // 					_tolerance,
        // 					_meanErrorFactor);

#ifdef HAVE_PKG_libprof

	ProfileEnd("interpBest");

#endif // HAVE_PKG_libprof

	return true;
	
      }
	
      //
      // this should never be reached
      //

      assert(false);

      return false;

    }

    double
    KrigingDataBase::interpolateSpecificModel(double            * value,
                                              double            * gradient,
                                              int               & model,
                                              const double      * point,
                                              std::vector<bool> & flags)
    {
       double errorEstimate;

       //
       // make sure there is enough space in flags 
       //

       assert(flags.size() >= NUMBER_FLAGS);

       //
       // initialize flags container
       //

       std::fill(flags.begin(),
                 flags.end(),
                 false);
       //
       // shortcuts to frequently accesses data
       //

       const int pointDimension = getPointDimension();
       const int valueDimension = getValueDimension();

       //
       // instatiate point object from point data
       //

       const ResponsePoint queryPoint(pointDimension,
                                      point);

       if (model == id_undefined) {
          cout << "Undefined model passed to adaptive sampler" << endl;
          exit(1);
       }

       uint128_t model_key = _ann.getKey(model);

       if (model_key == uint128_t_undefined) {

          cout << "Couldn't find model in kriging interpolation database " << endl;
          exit(1);

          //          flags[LOST_HINT_FLAG] = true;

       } else {
   
           InterpolationModelPtr hintKrigingModel = this->_modelDB->extract(model_key);
	  
          assert(hintKrigingModel->hasGradient() == true);

          flags[USED_HINT_FLAG] = true;

          errorEstimate = 0.;

          for (int iValue = 0; iValue < valueDimension; ++iValue) {

             const double iErrorEstimate = 
                compKrigingError(*hintKrigingModel, 
                                 queryPoint,
                                 iValue,
                                 _meanErrorFactor);
	  
             errorEstimate = std::max(errorEstimate, sqrt(iErrorEstimate));

             //
             // compute the value
             //

             const Value valueEstimate = 
                hintKrigingModel->interpolate(iValue,
                                          queryPoint);

             //
             // put the value of the function into value (valueEstimate
             // contains the value of the function follwed by the gradient;
             //

             value[iValue] = valueEstimate[0];

             //
             // store gradient data
             //
	  
             for (int i = 0; i < pointDimension; ++i)
                gradient[i*valueDimension + iValue] = valueEstimate[1 + i];

          }
       }

       return errorEstimate;
    }

    //
    // Insert the point-value pair into the database
    //

    void 
    KrigingDataBase::insert(int               & hint,
                            const double      * point,
                            const double      * value,
                            const double      * gradient,
                            std::vector<bool> & flags)
    {

      //
      // 
      // make sure there is enough space in flags 
      //

      assert(flags.size() >= NUMBER_FLAGS);

      //
      // initialize flags container
      //

      std::fill(flags.begin(),
		flags.end(),
		false);

#ifdef HAVE_PKG_libprof

      ProfileBegin("insert");

#endif // HAVE_PKG_libprof

      //
      // update number of point/value pairs
      //

      ++_numberPointValuePairs;

      //
      // shortcuts to frequently accessed data
      //

      const int pointDimension = getPointDimension();
      const int valueDimension = getValueDimension();

      //
      // check if the list of kriging models in non-empty
      //

      if (hint == id_undefined) {
	
	//
	// create and add new model
	//

         addNewModel(
                     _modelDB,
                     _ann,
                     _modelFactory,
                     hint, 
                     point,
                     value,
                     gradient,
                     pointDimension,
                     valueDimension);

	//
	// update number of kriging models
	//

	++_numberKrigingModels;

      } else {

	//
	// get a handle to the right kriging model
	//

         uint128_t model_key = _ann.getKey(hint);

        InterpolationModelPtr krigingModel = this->_modelDB->extract(model_key, (InterpolationModelFactoryPointer *)&this->_modelFactory);

	//
	// check the size of the model; if the next point would put
	// the model above _maxKrigingModelSize start a new model;
	// otherwise just add point/value pair
	//

	if (krigingModel->getNumberPoints() == _maxKrigingModelSize) {

           addNewModel(
                       _modelDB,
                       _ann,
                       _modelFactory,
                       hint, 
                       point,
                       value,
                       gradient,
                       pointDimension,
                       valueDimension);

	  //
	  // update number of kriging models
	  //

	  ++_numberKrigingModels;

	  //
	  // record event
	  //

	  flags[MODEL_SIZE_LIMIT_FLAG] = true;
 
	} else {

	  //
	  // copy value and gradient data
	  //
	  
	  std::vector<Value> pointValue;

	    
	  if (krigingModel->hasGradient() == true) 
	    pointValue = copyValueData(value,
				       gradient,
				       pointDimension,
				       valueDimension);
	  else
	    pointValue = copyValueData(value,
				       pointDimension,
				       valueDimension);
	  
	  //
	  // copy point data into a Point object
	  //
	  
	  const Point pointObject(pointDimension,
				  point);
	  
	  //
	  // insert point/value data into the model
	  //

	  const bool addPointSuccess = krigingModel->addPoint(pointObject,
							      pointValue);

	  if (addPointSuccess == true) {
	    
	    //
	    // compute the center of mass for the updated model
	    //
	    
	    const Point centerMass = getModelCenterMass(*krigingModel);

	    //
	    // remove old kriging model from the database
	    //

            _ann.remove(hint);

            _modelDB->erase(model_key);
	    
	    //
	    // insert updated kriging model into database
	    //
	    
	    const ResponsePoint centerMassRP(pointDimension,
					     &(centerMass[0]));
	    
            // Create a key string corresponding to the center of mass point
            // and insert it into the approximate nearest neighbor database

            std::vector<double> point_data;
            point_data.resize(centerMassRP.size());
            for (int i=0; i<centerMassRP.size(); ++i) {
               point_data[i] = centerMassRP[i];
            }

            uint128_t new_model_key = getKeyHash(centerMassRP);
            hint = _ann.insert(point_data, new_model_key);

            // Insert the interpolation model into the interpolation model database

            _modelDB->insert(new_model_key, krigingModel, (ResponsePoint *) &centerMassRP);

	  } else {

	    //
	    // point insertion failed-add new model
	    //

             addNewModel(
                         _modelDB,
                         _ann,
                         _modelFactory,
                         hint, 
                         point,
                         value,
                         gradient,
                         pointDimension,
                         valueDimension);

	    //
	    // update number of kriging models
	    //
	    
	    ++_numberKrigingModels;

	    //
	    // record event
	    //

	    flags[MODEL_INSERT_LIMIT_FLAG] = true;

	  }

	}

      }

#ifdef HAVE_PKG_libprof

      ProfileEnd("insert");

#endif // HAVE_PKG_libprof      
      
      return;

    }

    //
    // get the number of statistcs 
    //

    int
    KrigingDataBase::getNumberStatistics() const
    {

      return 2;

    }

    //
    // Write performance stats
    //

    void
    KrigingDataBase::getStatistics(double * stats,
                                   int      size) const
    {
      switch(size) {

      case 0:
	
	break;

      case 1:
	stats[0] = _numberKrigingModels;
	break;

      default:
	stats[0] = _numberKrigingModels;
	stats[1] = _numberPointValuePairs;
	break;

      }

      return;
    }

    //
    // Provide short descriptions of statistic data
    //

    std::vector<std::string>
    KrigingDataBase::getStatisticsNames() const
    {

      std::vector<std::string> names;
      
      names.push_back("Number of kriging models");
      names.push_back("Number of point/value pairs");

      return names;

    }

    //
    // output kriging model database stats
    //

    void
    KrigingDataBase::printDBStats(std::ostream & outputStream)
    {
      //
      // get number of statistics
      //

      const int numberStats = getNumberStatistics();
      
      //
      // get stats
      //

      std::vector<double> stats(numberStats);

      getStatistics(&(stats[0]),
		    numberStats);


      //
      // get descriptions of statistics
      //

      const std::vector<std::string> statStrings = getStatisticsNames();
      assert(static_cast<int>(statStrings.size()) == numberStats);
      
      //
      // output stats and description
      //

      for (int i = 0; i < numberStats; ++i)
	outputStream << statStrings[i] <<  " " << stats[i] << std::endl;

      //
      // output kriging model stats
      //

      //      outputKrigingModelStats(outputStream, _keyDB, _maxKrigingModelSize);
      
      //
      // output db stats
      //

      //      _keyDB.outputStats(outputStream);

      //
      // output positions of all kriging models
      //

      //      outputKrigingModelPositionData("kriging_model_centers.txt", _keyDB);

      return;
    }

    //
    // Swap out some objects to save memory
    //

    void
    KrigingDataBase::swapOutObjects() const
    {
#if 0
       if (typeid(_keyDB) == typeid(MTree)) {
          ((MTree&)_keyDB).writeObjects(KrigingModelChooser(_agingThreshold));
       }
       else {
          // Need to figure out how to generalize this to work
          // with databases other than MTree
       }
#else
       cout << "KrigingDataBase::swapOutObjects(): Figure out why this is being called" << endl;
       exit(1);
#endif

      return;

    }
}
