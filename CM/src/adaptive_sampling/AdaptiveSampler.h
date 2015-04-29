#ifndef ADAPTIVESAMPLER_INCLUDED
#define ADAPTIVESAMPLER_INCLUDED

#include <key_db/KrigingInterpolationKeyDB.h>
#include <vector>

#include "FineScale.h"

using namespace krigalg;
using namespace krigcpl;


class AdaptiveSampler
{
 public:

   AdaptiveSampler( const int                  pointDimension,
                    const int                  valueDimension,
                    const std::vector<double>& pointScaling,
                    const std::vector<double>& valueScaling,
                    const int                  maxKrigingModelSize,
                    const int                  maxNumberSearchModels,
                    const double               theta,
                    const double               meanErrorFactor,
                    const double               tolerance,
                    const double               maxQueryPointModelDistance );

   ~AdaptiveSampler();

   void sample( std::vector<double>&       value,
                const std::vector<double>& queryPoint,
                int&                       hint,
                std::vector<bool> &        flags,
                const FineScale&           fineScaleModel,
                double &                   error_estimate );

   void evaluateSpecificModel( std::vector<double>&       value,
                               const std::vector<double>& point,
                               int                        model,
                               const FineScale&           fineScaleModel,
                               double&                    error_estimate );

   void interpolateSpecificModel( std::vector<double>&       value,
                                  const std::vector<double>& point,
                                  int                        model,
                                  std::vector<bool>&         flags,
                                  double&                    error_estimate );

   void getModelInfo( int& numModels,
                      int& numPairs ) const;

   void printStatistics( std::ostream& outputStream );

   void printNewInterpolationStatistics( std::ostream& outputStream );

   int getNumSamples() const {return m_num_samples;}

   int getNumSuccessfulInterpolations() const {return m_num_successful_interpolations;}

   double getAveragePointNorm() const {return m_point_norm_sum / m_num_samples;}

   double getAverageValueNorm() const {return m_value_norm_sum / m_num_samples;}

   double getPointNormMax() const {return m_point_norm_max;}

   double getValueNormMax() const {return m_value_norm_max;}

   void setVerbose(const bool verbose) const {m_verbose = verbose;}

 private:

   void verifyInterpolationAccuracy( const std::vector<double>& point,
                                     const std::vector<double>& value,
                                     const FineScale&           fineScaleModel ) const;

   void printInterpolationFailure( const std::vector<double>& point,
                                   const std::vector<double>& value,
                                   const std::vector<double>& exact_value ) const;

   double pointL2Norm( const double* point ) const;

   double pointMaxNorm( const double* point ) const;

   double pointL2Norm( const std::vector<double>& point ) const;

   double pointMaxNorm( const std::vector<double>& point ) const;

   double valueL2Norm( const double* value ) const;

   double valueMaxNorm( const double* value ) const;

   double valueL2Norm( const std::vector<double>& value ) const;

   double valueMaxNorm( const std::vector<double>& value ) const;

   KrigingInterpolationKeyDB* m_interp;

   InterpolationModelDataBase m_modelDB;
   DB* m_keyDB;

   std::vector<double> m_pointScaling;
   std::vector<double> m_valueScaling;

   int m_pointDimension;
   int m_valueDimension;
   int m_valueAllocated;
   int m_maxKrigingModelSize;
   int m_maxNumberSearchModels;
   double m_theta;
   double m_meanErrorFactor;
   double m_tolerance;
   double m_maxQueryPointModelDistance;

   int m_printed_num_interp_models;
   int m_printed_num_pairs;
   int m_num_samples;
   int m_num_successful_interpolations;
   int m_num_fine_scale_evaluations;

   double m_point_norm_sum;
   double m_value_norm_sum;

   double m_point_norm_max;
   double m_value_norm_max;

   mutable bool m_verbose;
};

#endif
