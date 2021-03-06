#include "vpsc.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

extern "C"
{
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
         const int* diagnostics
         ); // populate all the required data structures 

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

         );  // call the VPSC fortran executable (with all data needed without any additional file accesses)
}

   void 
vpsc::printState()
{
   using std::cout;
   using std::endl;

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
vpsc::vpsc_init_class()
{
   char fnameIn[20]="vpsc_as_try.in";
   char str[1000];
   FILE *inFile=NULL;
   int nPhaseMax = 2;

   int iPhase, numPhases;

   inFile = fopen(fnameIn,"r");
   fscanf(inFile,"%d %[^\n]\n", &numPhases, str);
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
   nGrTot = 300;
   nH = 19;
   nTwinSysMax = 12;
   diagnostics = 1;
   enableUnitStress = 0;

   // allocate space for the materials
   // interface data
   hVecInit = (double*)malloc(nH * sizeof(double));
   
   // phase data
   nGrains     = (int*)malloc(nPhase * sizeof(int));
   nTwSysPhase = (int*)malloc(nPhase * sizeof(int));
   nStrPhase   = (int*)malloc(nPhase * sizeof(int));
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



   for (iPhase = 0; iPhase < nPhase; iPhase++) {
      fscanf(inFile, "%lf", &volFrac[iPhase]);
      printf("Phase %d has volume fraction %lf\n", iPhase, volFrac[iPhase]);
   }
   fscanf(inFile,"\n");

   fclose(inFile);
   
   printState();

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
         &diagnostics
         ); // populate all the required data structures 
   
   printState();

   printf("Finished with vpsc_init_\n");

}

double*
vpsc::tensorFunction(const double* in) const
{
   double* out;
   out = (double*)malloc(6*sizeof(double));

   double inFlat[6];
   double outFlat[6];

// copy tensor values to flat array
   for (int i = 0; i < 6; i++) { 
      inFlat[i] = in[i];
      printf("%f, %f\n", in[i], inFlat[i]);
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
      out[i] = outFlat[i];
   }
   
   return out;

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

/*
int main()
{
   vpsc test_state;

   test_state.vpsc_init_class();

   double inFlat[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

   inFlat[0] = 20.0*-1.0/3.0;
   inFlat[2] = 20.0*-1.0/3.0;
   inFlat[5] = 20.0* 2.0/3.0;

   test_state.tensorFunction(inFlat);
}
*/
