#ifndef _VPSC_
#define _VPSC_

#include "Plasticity.h"

class vpsc 
: public Plasticity
{
   public:

      vpsc(const double D_0,
            const double m,
            const double g,
            const double c_scaling)
         : m_D_0(D_0), m_m(m), m_g(g){vpsc_init_class(c_scaling);};

      ~vpsc(){;}

      virtual void vpsc_init_class(const double c_scaling);

      virtual Tensor2Sym tensorFunction( const Tensor2Sym& in ) const;

      //virtual void advance( const double delta_t, void* state);

      virtual void printState();

      virtual void getScalingsForSampling( vector<double>& input_scaling,
            vector<double>& output_scaling ) const;

      virtual void advance( const double delta_t, void* state);

      virtual void getState( void* state ) const;

      virtual size_t getStateSize() const {return 0;}

   private:

      double   m;

      double   g;

      // Taylor values for matching at the cutoff point
      //
      double m_D_0;      // Reference strain rate

      double m_m;        // Sensitivity

      double m_g;        // Hardness



      // phase_data_type
      int*     nGrains;
      double*  orients;
      double*  Cbin;
      double*  volFrac;
      double*  strengths;
      int*     iSymmCode;
      int*     nTwSysPhase;
      double*  twinOR;
      int*     iHVLB;
      int*     iHVUB;
      int*     nStrPhase;

      // phase_response_type
      double*  weights;
      double*  dw_dOdfDof;
      int*     nOdfDof;
      double*  twinVFRates;
      double*  reorientationRates;

      // interface_type
      int      iSCR;
      int      nPhase;
      int      nH;
      double   dRef;
      int      intPorosity;
      int      porosityIHVLB;
      int      porosityIHVUB;
      double*  hVecInit;
      bool     enableUnitStress;
      int      nTwinSysMax;
      int      nStrngthMax;
      int      nGrTot;

      int     diagnostics; // due to fortran - c compatibility
      double   s_threshold; 

};

#endif

