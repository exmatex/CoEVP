#ifndef _CONSTITUTIVE_
#define _CONSTITUTIVE_

#include "FineScale.h"
#include "tensor.h"

class AdaptiveSampler;

class Constitutive
{
 public:

   Constitutive()
      : m_sampler(NULL), m_finescale_verbose(false) {;}

   ~Constitutive();

   virtual void advance( const double delta_t ) = 0;

   virtual void setNewVelocityGradient( const Tensor2Gen& L_new ) = 0;

   virtual void setVolumeChange( const double volume_change ) = 0;

   virtual Tensor2Sym stress( const double compression,
                              const double e,
                              const double q ) const = 0;

   virtual Tensor2Sym stressDeviator() const = 0;

   virtual double pressure( const double compression,
                            const double internal_energy ) const = 0;

   virtual double soundSpeedSquared( const double reference_density,
                                     const double relativeVolume,
                                     const double energy ) const = 0;

   void enableAdaptiveSampling( const int                  pointDimension,
                                const int                  valueDimension,
                                const std::vector<double>& pointScaling,
                                const std::vector<double>& valueScaling,
                                const int                  maxKrigingModelSize,
                                const int                  maxNumberSearchModels,
                                const double               theta,
                                const double               meanErrorFactor,
                                const double               tolerance,
                                const double               maxQueryPointModelDistance );

   void sample( const FineScale&           fine_scale_model,
                const std::vector<double>& point,
                std::vector<double>&       value ) const;

   void evaluateSpecificModel( const int                  model,
                               const FineScale&           fine_scale_model,
                               const std::vector<double>& point,
                               std::vector<double>&       value ) const;

   bool adaptiveSamplingEnabled() const;

   void getModelInfo( int& numModels,
                      int& numPairs ) const;

   void printStats();

   void printNewInterpStats();

   int getNumSamples() const;

   int getNumSuccessfulInterpolations() const;

   double getAveragePointNorm() const;

   double getAverageValueNorm() const;

   double getPointNormMax() const;

   double getValueNormMax() const;

   mutable int m_hint;

   mutable double m_error_estimate;

 private:

   mutable bool m_finescale_verbose;

   AdaptiveSampler* m_sampler;

};

#endif

