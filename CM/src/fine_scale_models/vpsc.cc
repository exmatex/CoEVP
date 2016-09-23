#include "vpsc.h"
#include <stdlib.h>
#include <string.h>
extern "C"
{
   // populate all the required data structures 
   void vpsc_init_(
         // interface data
         const double* dRef, 
         const int* nPhase, 
         const int* nH, 
         const double* hVecInit,
         const int* porosityIHVLB,
         const int* porosityIHVUB,
         const int* nTwinSysMax, 
         const int* nStrngthMax, 
         const int* nGrTot,
         const int* intPorosity,
         // phase data
         const int* nGrains,
         const double* orients,
         const double* volFrac,
         const int* nStrPhase,
         const int* nTwSysPhase,
         const int* iHVLB, 
         const int* iHVUB,
         const double* strengths,
         // phase response data
         const double* weights,
         const double* twinVFRates,
         const double* reorientationRates,
         const double* c_scaling,
         const int* diagnostics
            ); 
   // call the VPSC fortran executable (with all data needed without any additional file accesses)
   void vpsc_run_(
         double* in, 
         double* out,
         const double* m, 
         const double* g, 
         const double* dRef, 
         const int* nPhase, 
         const int* nH, 
         const double* hVecInit,
         const int* porosityIHVLB,
         const int* porosityIHVUB,
         const int* nTwinSysMax, 
         const int* nStrngthMax, 
         const int* nGrTot,
         const int* intPorosity,
         // phase data
         const int* nGrains,             
         const double* orients,       
         const double* volFrac,       
         // shapes,        
         const int* nStrPhase,     
         const int* nTwSysPhase,    
         const int* iHVLB, 
         const int* iHVUB,        
         const double* strengths,     
         // phase response data
         const double* weights,             
         const double* twinVFRates,         
         const double* reorientationRates,  
         const int* diagnostics          

            );  
}

   void 
vpsc::printState()
{
   int i, j;
   int iPhase, iGlobal, iWeights;

   cout << "m = " << m << endl;
   cout << "g = " << g << endl;

   cout << "Interface variables" << endl;

   cout << " iSCR = " << iSCR << endl;
   cout << " nPhase = " << nPhase << endl;
   cout << " nH = " << nH << endl;
   cout << " dRef = " << dRef << endl;
   cout << " intPorosity = " << intPorosity << endl;
   cout << " porosityIHVLB = " << porosityIHVLB << endl;
   cout << " porosityIHVUB = " << porosityIHVUB << endl;
   cout << " enableUnitStress = " << enableUnitStress << endl;
   cout << " nTwinSysMax = " << nTwinSysMax << endl;
   cout << " nStrngthMax = " << nStrngthMax << endl;
   cout << " nGrTot = " << nGrTot << endl;
   cout << " hVecInit: "<< endl;
   for (i = 0; i < nH; i++) {
      cout << hVecInit[i] << endl;
   }

   iGlobal = 0;
   for (iPhase = 0; iPhase < nPhase; iPhase++) {
      cout << "Phase Data " << iPhase << endl;
      cout << "   nGrains = " << nGrains[iPhase] << endl;
      cout << "   volFrac = " << volFrac[iPhase] << endl;
      cout << "   iSymmCode = " << iSymmCode[iPhase] << endl;
      cout << "   nTwinSystems = " << nTwSysPhase[iPhase] << endl;
      cout << "   iHVLB = " << iHVLB[iPhase] << endl;
      cout << "   iHVUB = " << iHVUB[iPhase] << endl;
      cout << "   nStrengths = " << nStrPhase[iPhase] << endl;
   }
   iGlobal = 0;
   cout << "   strengths: " << endl << endl;
   for (iPhase = 0; iPhase < nPhase; iPhase++) {
      for (i = 0; i < nStrPhase[iPhase]; i++) {
         printf("%10.4f\n", 
               strengths[iGlobal]);
         iGlobal++ ;
      }
   }
   iGlobal = 0;
   cout << "   orientations, weights: " << endl << endl;
   for (iPhase = 0; iPhase < nPhase; iPhase++) {
      for (i = 0; i < nGrains[iPhase]; i++) {
         printf("%10.4f, %10.4f, %10.4f, %10.4f\n", 
               orients[3*iGlobal+0], orients[3*iGlobal+1], orients[3*iGlobal+2], weights[iGlobal]);
         iGlobal++ ;
      }
   }
   fflush(stdout);
}

   void
vpsc::vpsc_init_class(const double c_scaling)
{
   diagnostics = 0;
   char str[1000];
   char fnameIn[255];
   char dataDir[255];
   char threshold[255];

   if(std::getenv("VPSC_INPUT_PATH")==NULL)
   {
	   strcpy(fnameIn, "../../CoEVP/CM/src/fine_scale_models/tantalum/vpsc_as_try.in");
   }
   else
   {
	  strcpy(fnameIn,std::getenv("VPSC_INPUT_PATH"));
	  strcat(fnameIn,"vpsc_as_try.in");
   }
   //printf("Opening %s\n", fnameIn);

   FILE *inFile=NULL;
   int nPhaseMax = 2;

   int iPhase, numPhases;

   inFile = fopen(fnameIn,"r");
   fscanf(inFile,"%d %[^\n]\n", &numPhases, str);
   if (diagnostics == 1)
      printf("numPhases = %d\n", numPhases);
   if (numPhases > nPhaseMax) {
      printf("We don't handle more than %d phases\n", nPhaseMax);
      exit;
   }
   // parameters from the fortran code
   int eulerdim_p = 3;
   int maxStrengthsPerPhase = 10;
   // kluge for now
   nPhase = numPhases;
   nGrTot = 500;
   nH = 19;
   nTwinSysMax = 12;

   // allocate space for the materials
   // interface data
   hVecInit = (double*)malloc(nH * sizeof(double));

   // phase data
   nGrains     = (int*)malloc(nPhase * sizeof(int));
   nStrPhase   = (int*)malloc(nPhase * sizeof(int));
   nTwSysPhase = (int*)malloc(nPhase * sizeof(int));
   iHVLB       = (int*)malloc(nPhase * sizeof(int));
   iHVUB       = (int*)malloc(nPhase * sizeof(int));
   iSymmCode   = (int*)malloc(nPhase * sizeof(int));
   volFrac     = (double*)malloc(nPhase * sizeof(double));
   orients     = (double*)malloc(eulerdim_p*nGrTot * sizeof(double));
   strengths   = (double*)malloc(nPhase*maxStrengthsPerPhase* sizeof(double));
   // phase response data
   weights     = (double*)malloc(nGrTot * sizeof(double));
   twinVFRates = (double*)malloc(nTwinSysMax*nGrTot * sizeof(double));
   reorientationRates = (double*)malloc(eulerdim_p*nGrTot * sizeof(double));
   // stress threshold for bypass 
   if(std::getenv("VPSC_S_THRESHOLD") == NULL) { 
     s_threshold = 1.0e-4; 
   } else { 
     strcpy(threshold,std::getenv("VPSC_S_THRESHOLD"));
     s_threshold  = strtod(threshold, NULL); 
     if(diagnostics) 
       std::cout << "VPSC Init: set stress threshold to: " << threshold
                 << " double value is : " << s_threshold << '\n';
   }
   for (iPhase = 0; iPhase < nPhase; iPhase++) {
      fscanf(inFile, "%lf", &volFrac[iPhase]);
   if (diagnostics == 1)
      printf("Phase %d has volume fraction %lf\n", iPhase, volFrac[iPhase]);
   }
   fscanf(inFile,"\n");

   fclose(inFile);

   //printState();

   if (diagnostics == 1)
   printf("Calling vpsc_init_\n");

   // initialize the values using the fortran routines
   vpsc_init_(
         // interface data
         &dRef, 
         &nPhase, 
         &nH, 
         hVecInit,
         &porosityIHVLB,
         &porosityIHVUB,
         &nTwinSysMax, 
         &nStrngthMax, 
         &nGrTot,
         &intPorosity,
         // phase data
         nGrains,
         orients,
         volFrac,
         nStrPhase,
         nTwSysPhase,
         iHVLB, 
         iHVUB,
         strengths,
         // phase response data
         weights,
         twinVFRates,
         reorientationRates,
         &c_scaling,
         &diagnostics
            ); // populate all the required data structures 

   //printState();

   if (diagnostics == 1)
   printf("Finished with vpsc_init_\n");

}

Tensor2Sym
vpsc::tensorFunction(const Tensor2Sym& in) const
{
   Tensor2Sym out;
   Tensor2Sym in_dev;

   in_dev=dev(in);

   double normThreshold =  s_threshold;
   double inFlat[6];
   double outFlat[6];

   double normIn = norm(in_dev);

   // bypass vpsc call if stress is too small
   if (normIn > normThreshold) 
   {
      // copy tensor values to flat array
      for (int i = 0; i < 6; i++) { 
         inFlat[i] = in_dev.a[i];
         if (diagnostics > 0)
         printf("%f, %f\n", in_dev.a[i], inFlat[i]);
      }

      vpsc_run_(
            inFlat, 
            outFlat,
            // interface variables
            &m, &g, &dRef, 
            &nPhase, 
            &nH, 
            hVecInit,
            &porosityIHVLB,
            &porosityIHVUB,
            &nTwinSysMax, 
            &nStrngthMax, 
            &nGrTot,
            &intPorosity,
            // phase data
            nGrains,             
            orients,       
            volFrac,       
            // shapes,        
            nStrPhase,     
            nTwSysPhase,    
            iHVLB, 
            iHVUB,        
            strengths,     
            // phase response data
            weights,             
            twinVFRates,         
            reorientationRates,  
            &diagnostics          
               );

      // copy back to output tensor
      for (int i = 0; i < 6; i++) { 
         out.a[i] = outFlat[i];
      }
   } else {
      printf("VPSC call bypass\n");
      for (int i = 0; i < 6; i++) { 
         //in.a[i] = 0.0; // zero the input also for consistency
         //out.a[i] = 0.0;
         out.a[i] = 1.0e2*in_dev.a[i];
         printf("%g, %g\n", in_dev.a[i], out.a[i]);
      }
   }
   /*
   for (int i = 0; i < 6; i++) { 
       printf("%g, %g\n", in.a[i], out.a[i]);
   }
   */

   return out;

}

void
vpsc::getScalingsForSampling( vector<double>& input_scaling,
      vector<double>& output_scaling ) const
{
   assert(input_scaling.size() == m_pointDimension);
   assert(output_scaling.size() == m_valueDimension);

   /*
      if ( m == 1. ) {
      for (int i=0; i<m_pointDimension; ++i) {
      input_scaling[i] = 1.e-1;
      }
      for (int i=0; i<m_valueDimension; ++i) {
      output_scaling[i] = 1.e3;
      }
      }
      else if ( m == 1./2. ) {
      for (int i=0; i<m_pointDimension; ++i) {
   //         input_scaling[i] = 1.e-2;
   input_scaling[i] = 1.e-4;
   }
   for (int i=0; i<m_valueDimension; ++i) {
   output_scaling[i] = 1.e1;
   }
   }
   else {
   */
   for (int i=0; i<m_pointDimension; ++i) {
      input_scaling[i] = 1.0;
   }
   for (int i=0; i<m_valueDimension; ++i) {
      output_scaling[i] = 1.0;
   }
   //}
}


   void
vpsc::advance( const double delta_t, void* state)
{
   // This doesn't actually do anything but copy three doubles
   // back and forth, but is here as an example of how a more
   // sophisticated plasticity model might work

   // Get the state
   double* local_state = (double*)state;

   // Advance the fine-scale model (trivially)

   // Update the state
   local_state = (double*)state;

}


void
vpsc::getState( void* state ) const
{
   double* local_state = (double*)state;
}

