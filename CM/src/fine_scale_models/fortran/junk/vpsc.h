#ifndef _VPSC_
#define _VPSC_


class vpsc 
{
   public:

      vpsc() {};

      ~vpsc(){;}

      //{ vpsc_init_class();  };

      virtual void vpsc_init_class();

      virtual double* tensorFunction( const double* in ) const;

      //virtual void advance( const double delta_t, void* state);

      virtual void printState();

      virtual void advance( const double delta_t, void* state);

      virtual void getState( void* state ) const;

   private:

      double   m;
      double   g;


      // phase_data_type
      int*     nGrains;
      int*     nStrPhase;
      int*     nTwSysPhase;
      int*     iHVLB;
      int*     iHVUB;
      int*     iSymmCode;
      double*  orients;
      double*  Cbin;
      double*  volFrac;
      double*  strengths;
      double*  twinOR;

      // phase_response_type
      double*  weights;
      //double*  dw_dOdfDof;
      //int*     nOdfDof;
      double*  twinVFRates;
      double*  reorientationRates;

      // interface_type
      int      iSCR;
      int      nPhase;
      int      nH;
      int      intPorosity;
      int      porosityIHVLB;
      int      porosityIHVUB;
      int      nTwinSysMax;
      int      nStrngthMax;
      int      nGrTot;
      double   dRef;
      double*  hVecInit;
      // due to fortran - c compatibility bools are cast to int
      int     enableUnitStress;
      int     diagnostics; 
};

#endif

