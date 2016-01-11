#ifndef _CONSTITUTIVE_
#define _CONSTITUTIVE_

#include <cstring>
#include "ConstitutiveGlobal.h"
#include "ApproxNearestNeighbors.h"
#include "AdaptiveSampler.h"
#include "FineScale.h"
#include "tensor.h"
#include "ModelDatabase.h"

struct ConstitutiveData
{
   Tensor2Sym  sigma_prime;
   int         num_models;
   int         num_point_value_pairs;
   int         num_Newton_iters;
};

class Constitutive
{
 public:

   Constitutive( ConstitutiveGlobal& global)
      : m_global(&global), m_sampler(NULL), m_finescale_verbose(false) {;}

   ~Constitutive();

   //   virtual void advance( const double delta_t ) = 0;
   virtual ConstitutiveData advance( const double delta_t, const Tensor2Gen& L_new, const double, void* state ) = 0;

   virtual Tensor2Sym stress( const double compression,
                              const double e,
                              const double q ) const = 0;

   virtual Tensor2Sym stressDeviator() const = 0;

   virtual double pressure( const double compression,
                            const double internal_energy ) const = 0;

   virtual double soundSpeedSquared( const double reference_density,
                                     const double relativeVolume,
                                     const double energy ) const = 0;

   virtual size_t getStateSize() const = 0;

   virtual void getState( void* buffer ) const = 0;

   virtual void* setState( void* buffer ) = 0;

   void enableAdaptiveSampling( const int                  pointDimension,
                                const int                  valueDimension,
                                const std::vector<double>& pointScaling,
                                const std::vector<double>& valueScaling,
                                const int                  maxKrigingModelSize,
                                const int                  maxNumberSearchModels,
                                const double               theta,
                                const double               meanErrorFactor,
                                const double               tolerance,
                                const double               maxQueryPointModelDistance,
                                ApproxNearestNeighbors*    ann,
                                ModelDatabase*             modelDB);

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

 protected:
   
   template <class T>
   inline void pushObjectToBuffer( void**     buffer,
                                   const T&   data ) const
   {
      size_t object_size = sizeof(T);
      memcpy(*buffer, &data, object_size);
      *buffer = ((char*)(*buffer)) + object_size;
   };

   template <class T>
   inline void popObjectFromBuffer( void** buffer,
                                    T&     data )
   {
      size_t object_size = sizeof(T);
      memcpy(&data, *buffer, object_size);
      *buffer = ((char*)(*buffer)) + object_size;
   };

 public:

   mutable bool m_finescale_verbose;

   AdaptiveSampler* m_sampler;

   ConstitutiveGlobal* m_global;

};


#endif

