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

#include "Constitutive.h"
#include "AdaptiveSampler.h"


Constitutive::~Constitutive()
{
   if ( adaptiveSamplingEnabled() ) {
      delete m_sampler;
      m_sampler = NULL;
   }
}


bool
Constitutive::adaptiveSamplingEnabled() const
{
   return m_sampler != NULL;
}


void
Constitutive::enableAdaptiveSampling( const int                  pointDimension,
                                      const int                  valueDimension,
                                      const std::vector<double>& pointScaling,
                                      const std::vector<double>& valueScaling,
                                      const int                  maxKrigingModelSize,
                                      const int                  maxNumberSearchModels,
                                      const double               theta,
                                      const double               meanErrorFactor,
                                      const double               tolerance,
                                      const double               maxQueryPointModelDistance )
{
   m_sampler = new AdaptiveSampler( pointDimension,
                                    valueDimension,
                                    pointScaling,
                                    valueScaling,
                                    maxKrigingModelSize,
                                    maxNumberSearchModels,
                                    theta,
                                    meanErrorFactor,
                                    tolerance,
                                    maxQueryPointModelDistance );

   // This variable remembers the index of the most recently used interpolation model
   m_hint = -1;

   // This variable remembers the error estimate of the most recently used interpolation
   m_error_estimate = 0.;
}


void
Constitutive::sample( const FineScale&           fine_scale_model,
                      const std::vector<double>& point,
                      std::vector<double>&       value ) const
{
   std::vector<bool> interpolateFlags(InterpolationDataBase::NUMBER_FLAGS);

   m_sampler->setVerbose(true);

   m_sampler->sample( value,
                      point,
                      m_hint,
                      interpolateFlags,
                      fine_scale_model,
                      m_error_estimate );
}


void
Constitutive::evaluateSpecificModel( const int                  model,
                                     const FineScale&           fine_scale_model,
                                     const std::vector<double>& point,
                                     std::vector<double>&       value ) const
{
   m_sampler->evaluateSpecificModel( value,
                                     point,
                                     model,
                                     fine_scale_model,
                                     m_error_estimate );
}


void
Constitutive::getModelInfo( int& numModels, int& numPairs ) const
{
   if ( adaptiveSamplingEnabled() ) {
      m_sampler->getModelInfo(numModels, numPairs);
   }
   else {
      numModels = numPairs = 0;
   }
}


void
Constitutive::printStats()
{
   if ( adaptiveSamplingEnabled() ) {
      m_sampler->printStatistics(std::cout);
   }
}


void
Constitutive::printNewInterpStats()
{
   if ( adaptiveSamplingEnabled() ) {
      m_sampler->printNewInterpolationStatistics(std::cout);
   }
}


int
Constitutive::getNumSamples() const
{
   int samples = 0;

   if ( adaptiveSamplingEnabled() ) {
      samples = m_sampler->getNumSamples();
   }
   else {
      cout << "Constitutive::getNumSamples(): Adaptive sampling is not enabled!" << endl;
   }

   return samples;
}


int
Constitutive::getNumSuccessfulInterpolations() const
{
   int interps = 0;

   if ( adaptiveSamplingEnabled() ) {
      interps = m_sampler->getNumSuccessfulInterpolations();
   }
   else {
      cout << "Constitutive::getNumSuccessfulInterpolations(): Adaptive sampling is not enabled!" << endl;
   }

   return interps;
}


double
Constitutive::getAveragePointNorm() const
{
   double average = 0.;

   if ( adaptiveSamplingEnabled() ) {
      average = m_sampler->getAveragePointNorm();
   }

   return average;
}


double
Constitutive::getAverageValueNorm() const
{
   double average = 0.;

   if ( adaptiveSamplingEnabled() ) {
      average = m_sampler->getAverageValueNorm();
   }

   return average;
}


double
Constitutive::getPointNormMax() const
{
   double max = 0.;

   if ( adaptiveSamplingEnabled() ) {
      max = m_sampler->getPointNormMax();
   }

   return max;
}


double
Constitutive::getValueNormMax() const
{
   double max = 0.;

   if ( adaptiveSamplingEnabled() ) {
      max = m_sampler->getValueNormMax();
   }

   return max;
}

