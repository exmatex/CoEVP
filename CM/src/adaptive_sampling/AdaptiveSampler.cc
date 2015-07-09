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

#include "AdaptiveSampler.h"

#include <kriging/LinearDerivativeRegressionModel.h>
#include <kriging/GaussianDerivativeCorrelationModel.h>
#include <kriging/MultivariateDerivativeKrigingModelFactory.h>

#include <mtreedb/MTree.h>

#include "interpolation_database/kriging_db/DBKrigingModelObjectFactory.h"

AdaptiveSampler::AdaptiveSampler( const int                  pointDimension,
                                  const int                  valueDimension,
                                  const std::vector<double>& pointScaling,
                                  const std::vector<double>& valueScaling,
                                  const int                  maxKrigingModelSize,
                                  const int                  maxNumberSearchModels,
                                  const double               theta,
                                  const double               meanErrorFactor,
                                  const double               tolerance,
                                  const double               maxQueryPointModelDistance )
   : m_pointDimension(pointDimension),
     m_valueDimension(valueDimension),
     m_pointScaling(pointScaling),
     m_valueScaling(valueScaling),
     m_valueAllocated((pointDimension+1)*valueDimension),
     m_maxKrigingModelSize(maxKrigingModelSize),
     m_maxNumberSearchModels(maxNumberSearchModels),
     m_theta(theta),
     m_meanErrorFactor(meanErrorFactor),
     m_tolerance(tolerance),
     m_maxQueryPointModelDistance(maxQueryPointModelDistance),
     m_printed_num_interp_models(0),
     m_printed_num_pairs(0),
     m_num_samples(0),
     m_num_successful_interpolations(0),
     m_num_fine_scale_evaluations(0),
     m_point_norm_sum(0.),
     m_value_norm_sum(0.),
     m_point_norm_max(0.),
     m_value_norm_max(0.),
     m_verbose(false)
{

   DerivativeRegressionModelPointer 
      regressionModel(new LinearDerivativeRegressionModel);

   DerivativeCorrelationModelPointer
      correlationModel(new GaussianDerivativeCorrelationModel(std::vector<double>(1, m_theta)));

   MultivariateDerivativeKrigingModelFactoryPointer 
      modelFactory(new MultivariateDerivativeKrigingModelFactory(regressionModel,
                                                                 correlationModel));

   // Construct the key database

   m_keyDB = (DB*)(new MTree("kriging_model_database", &(std::cout), false));

   std::string mtreeDirectoryName = ".";

   m_keyDB->initializeCreate(mtreeDirectoryName + "/" 
                          "kriging_model_database",
                          "krigcpl",
                          *(new DBKrigingModelObjectFactory<InterpolationModel>));
      
   ((MTree*)m_keyDB)->setMaxNodeEntries(12);
     
   bool db_from_file = false;  // FIX THIS (input from somewhere)

   if ( db_from_file ) {
      cout << "AdaptiveSampler.cc: Database read option not yet implemented" << endl;
   }
   else {

      m_interp = new KrigingInterpolationKeyDB( m_pointDimension,
                                                m_valueDimension,
                                                modelFactory,
                                                *m_keyDB,
                                                m_modelDB,
                                                m_maxKrigingModelSize,
                                                m_maxNumberSearchModels,
                                                true,
                                                m_meanErrorFactor,
                                                m_tolerance,
                                                m_maxQueryPointModelDistance,
                                                600000000 );
   }

   // Add the gradient scaling from the point and value scaling
   m_valueScaling.resize(m_valueAllocated);
   for ( int i=m_valueDimension; i<m_valueAllocated; ++i ) {

      const int pointId = (i - m_valueDimension)/m_valueDimension;
      const int valueId = (i - m_valueDimension) - pointId*m_valueDimension;
	
      assert(pointId >=0 && pointId < m_pointDimension);
      assert(valueId >=0 && valueId < m_valueDimension);

      m_valueScaling[i] = m_valueScaling[valueId] / m_pointScaling[pointId];
   }
}


AdaptiveSampler::~AdaptiveSampler()
{
   m_pointScaling.resize(0);
   m_valueScaling.resize(0);
   
   delete m_interp;
   delete m_keyDB;
}


void
AdaptiveSampler::sample( std::vector<double>&       value,
                         const std::vector<double>& point,
                         int&                       hint,
                         std::vector<bool>&         flags,
                         const FineScale&           fineScaleModel,
                         double&                    error_estimate)
{
   // The interpolation database takes pointers to doubles rather
   // than stl vectors, so we need to make temporaries for the
   // query point and returned value.
   
   int point_length = point.size();

   double* local_point = new double[point_length];
   for (int i=0; i<point_length; ++i) {
      local_point[i] = point[i] / m_pointScaling[i];
   }

   double local_point_norm = pointL2Norm(local_point);
   if ( local_point_norm > m_point_norm_max ) {
      m_point_norm_max = local_point_norm;
   }
   m_point_norm_sum += local_point_norm;

   int value_length = value.size();
   double* local_value = new double[value_length];

   bool interpolationSuccess = 
      m_interp->interpolate(local_value,
                            hint,
                            local_point,
                            flags,
                            error_estimate);

   if (interpolationSuccess == false) {

      fineScaleModel.evaluate(point, value);

      if (m_verbose) {
#if 0
         std::cout << "   Adding key : (" << local_point[0];
         for (int i=1; i<m_pointDimension; ++i) {
            std::cout << ", " << local_point[i];
         }
         std::cout << ")";

         std::cout << ", value : (" << value[0] / m_valueScaling[0];
         for (int i=1; i<m_valueDimension; ++i) {
            std::cout << ", " << value[i] / m_valueScaling[i];
         }
         std::cout << ")" << std::endl;
#endif
         //         cout << "Interpolation failed: Adding key" << endl;
      }

      for (int i=0; i<value_length; ++i) {
         local_value[i] = value[i] / m_valueScaling[i];
      }

      m_interp->insert(hint,
                       local_point,
                       local_value,
                       &(local_value[m_valueDimension]),
                       flags);

      m_num_fine_scale_evaluations++;

      hint = -1;
      error_estimate = 0.;

   }
   else {

      if (m_verbose) {
         //         cout << "Interpolation succeeded" << endl;
      }

      interpolationSuccess = 
      m_interp->interpolate(local_value,
                            local_value + point_length,
                            hint,
                            local_point,
                            flags,
                            error_estimate);

      for (int i=0; i<value_length; ++i) {
         value[i] = local_value[i] * m_valueScaling[i];
      }

#if 0
      // For debugging: override the interpolant with the true fine-scale
      // model value and/or derivative

      std::vector<double> tmp_point;
      tmp_point.resize(m_pointDimension);
      for (int i=0; i<point_length; ++i) tmp_point[i] = point[i];

      std::vector<double> tmp_value;
      tmp_value.resize(value_length);

      fineScaleModel.evaluate(tmp_point, tmp_value);
      for (int i=0; i<value_length; ++i) value[i] = tmp_value[i];
#endif

      //      verifyInterpolationAccuracy(point, value, fineScaleModel);

      m_num_successful_interpolations++;
   }

   m_num_samples++;

   assert(m_num_samples = m_num_fine_scale_evaluations + m_num_successful_interpolations);

   double local_value_norm = valueL2Norm(local_value);
   if ( local_value_norm > m_value_norm_max ) {
      m_value_norm_max = local_value_norm;
   }
   m_value_norm_sum += local_value_norm;

   delete [] local_value;
   delete [] local_point;
}


void
AdaptiveSampler::evaluateSpecificModel( std::vector<double>&       value,
                                        const std::vector<double>& point,
                                        int                        model,
                                        const FineScale&           fineScaleModel,
                                        double&                    error_estimate )
{
   if ( model < 0 ) {
      fineScaleModel.evaluate(point, value);

      error_estimate = 0.;

      m_num_fine_scale_evaluations++;
   }
   else {

      // The interpolation database takes pointers to doubles rather
      // than stl vectors, so we need to make temporaries for the
      // query point and returned value.

      int point_length = point.size();

      double* local_point = new double[point_length];
      for (int i=0; i<point_length; ++i) {
         local_point[i] = point[i] / m_pointScaling[i];
      }

      int value_length = value.size();
      double* local_value = new double[value_length];

      std::vector<bool> interpolateFlags(InterpolationDataBase::NUMBER_FLAGS);

      error_estimate = m_interp->interpolateSpecificModel(local_value,
                                                      local_value + point_length,
                                                      model,
                                                      local_point,
                                                      interpolateFlags);

      for (int i=0; i<value_length; ++i) {
         value[i] = local_value[i] * m_valueScaling[i];
      }

      delete [] local_value;
      delete [] local_point;

      // Since a specific model is being used, the interpolation is deemed successful
      // by definition.
      m_num_samples++;
      m_num_successful_interpolations++;
   }
}


void
AdaptiveSampler::verifyInterpolationAccuracy( const std::vector<double>& point,
                                              const std::vector<double>& value,
                                              const FineScale&           fineScaleModel ) const
{
   // This function performs an independent test of the accuracy of the
   // kriging interpolation.  Only the function value is checked; not the
   // derivatives.  Since it makes a fine scale model call to compute
   // the true error, which is exactly what we hoped to avoid by doing
   // adaptive sampling, this function should obviously only be used
   // for verification or debugging purposes, e.g., if there is some
   // question about the effectiveness of the kriging error estimate.

   std::vector<double> exact_value(m_valueAllocated);
   fineScaleModel.evaluate(point, exact_value);

   double exact_norm = valueMaxNorm(exact_value);
   double eps = 1.e-11;

   if (exact_norm > eps) {

      std::vector<double> error(m_valueDimension);
      for (int i=0; i<m_valueDimension; ++i) {
         error[i] = value[i] - exact_value[i];
      }
      double error_norm = valueMaxNorm(error) / exact_norm;

      if (error_norm > m_tolerance) {
         cout << "error norm = " << error_norm << endl;
         //         printInterpolationFailure(point, value, exact_value);
         //         exit(1);
      }
   }
   else {  // Perform an absolute comparison if exact norm is small

      double value_norm = valueMaxNorm(value);

      //      if (value_norm > (1. + m_tolerance)*eps) {
      if (value_norm > m_tolerance) {
         cout << "value norm = " << value_norm << endl;
         //         printInterpolationFailure(point, value, exact_value);
         //         exit(1);
      }
   }
}


void
AdaptiveSampler::getModelInfo( int& numModels,
                               int& numPairs ) const
{
   double stats[2];
   m_interp->getStatistics(stats, 2);

   numModels = stats[0];
   numPairs = stats[1];
}


void
AdaptiveSampler::printStatistics( std::ostream & outputStream )
{
   m_interp->printDBStats(outputStream);
}



void
AdaptiveSampler::printNewInterpolationStatistics( std::ostream & outputStream )
{
   double stats[2];
   m_interp->getStatistics(stats, 2);

   bool print_stats = false;

   if (stats[0] > m_printed_num_interp_models) {
      m_printed_num_interp_models = stats[0];
      print_stats = true;
   }

   if (stats[1] > m_printed_num_pairs) {
      m_printed_num_pairs = stats[1];
      print_stats = true;
   }

   if (print_stats) {
      //      m_interp->printDBStats(outputStream);
      outputStream << "   # samples = " << m_num_samples
                   << ", # fine scale evaluations = " << m_num_fine_scale_evaluations
                   << ", successful interpolation ratio = "
                   << (double)m_num_successful_interpolations / (double)m_num_samples << endl;
   }
}


void
AdaptiveSampler::printInterpolationFailure( const std::vector<double>& point,
                                            const std::vector<double>& value,
                                            const std::vector<double>& exact_value ) const
{
   std::cout << "Interpolation failure at (" << point[0];
   for (int i=1; i<m_pointDimension; ++i) {
      std::cout << ", " << point[i];
   }
   std::cout << ")";

   std::cout << ", value = (" << value[0];
   for (int i=1; i<m_valueDimension; ++i) {
      std::cout << ", " << value[i];
   }
   std::cout << ")";

   std::cout << ", exact = (" << exact_value[0];
   for (int i=1; i<m_valueDimension; ++i) {
      std::cout << ", " << exact_value[i];
   }
   std::cout << ")" << endl;
}


inline double
AdaptiveSampler::pointL2Norm( const double* point ) const
{
   double norm = 0.;
   for (int i=0; i<m_pointDimension; ++i) {
      norm += point[i] * point[i];
   }

   return sqrt(norm);
}


inline double
AdaptiveSampler::pointMaxNorm( const double* point ) const
{
   double norm = 0.;
   for (int i=0; i<m_pointDimension; ++i) {
      double local_norm = fabs(point[i]);
      if (local_norm > norm) norm = local_norm;
   }

   return norm;
}


inline double
AdaptiveSampler::pointL2Norm( const std::vector<double>& point ) const
{
   assert(m_pointDimension <= point.size());

   double norm = 0.;
   for (int i=0; i<m_pointDimension; ++i) {
      norm += point[i] * point[i];
   }

   return sqrt(norm);
}


inline double
AdaptiveSampler::pointMaxNorm( const std::vector<double>& point ) const
{
   assert(m_pointDimension <= point.size());

   double norm = 0.;
   for (int i=0; i<m_pointDimension; ++i) {
      double local_norm = fabs(point[i]);
      if (local_norm > norm) norm = local_norm;
   }

   return norm;
}


inline double
AdaptiveSampler::valueL2Norm( const double* value ) const
{
   double norm = 0.;
   for (int i=0; i<m_valueDimension; ++i) {
      norm += value[i] * value[i];
   }

   return sqrt(norm);
}


inline double
AdaptiveSampler::valueMaxNorm( const double* value ) const
{
   double norm = 0.;
   for (int i=0; i<m_valueDimension; ++i) {
      double local_norm = fabs(value[i]);
      if (local_norm > norm) norm = local_norm;
   }

   return norm;
}


inline double
AdaptiveSampler::valueL2Norm( const std::vector<double>& value ) const
{
   assert(m_valueDimension <= value.size());

   double norm = 0.;
   for (int i=0; i<m_valueDimension; ++i) {
      norm += value[i] * value[i];
   }

   return sqrt(norm);
}


inline double
AdaptiveSampler::valueMaxNorm( const std::vector<double>& value ) const
{
   assert(m_valueDimension <= value.size());

   double norm = 0.;
   for (int i=0; i<m_valueDimension; ++i) {
      double local_norm = fabs(value[i]);
      if (local_norm > norm) norm = local_norm;
   }

   return norm;
}

