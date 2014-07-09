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
// File:	MemoryUtilities.cc
// Package:	toolbox
// 
// 
// 
//

#include "toolbox/base/MemoryUtilities.h"

#include "toolbox/parallel/MPI.h"
#include "toolbox/base/Utilities.h"
#include "toolbox/base/MathUtilities.h"
#include "toolbox/stream/IOStream.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_TAU
#include <malloc.h>
#endif

#ifdef HAVE_TAU
#if (PROFILING_ON || TRACING_ON)
#include <Profile/Profiler.h>
/* Register an "event" with Tau to track memory usage. */
TAU_PROFILE_STMT(TauUserEvent ue("memory use"))
#endif
#endif


namespace toolbox {

double MemoryUtilities::s_max_memory = 0.;

/*
*************************************************************************
*                                                                       *
* Prints memory usage to specified output stream.  Each time this       *
* method is called, it prints in the format:                            *
*                                                                       *
*    253.0MB (265334688) in 615 allocs, 253.9MB reserved (871952 unused)*
*                                                                       *
* where                                                                 *
*                                                                       *
*    253.0MB is how many megabytes your current allocation has malloced.* 
*    2653346688 is the precise size (in bytes) of your current alloc.   *
*    615 is the number of items allocated with malloc.                  *
*    253.9MB is the current memory reserved by the system for mallocs.  *
*    871952 is the bytes currently not used in this reserved memory.    *
*                                                                       *
*************************************************************************
*/
void MemoryUtilities::printMemoryInfo(ostream& os) 
{
#ifdef HAVE_TAU
   /*
    * NOTE: This was taken directly from John Gyllenhal...
    */
   
   /* Get malloc info structure */
   struct mallinfo my_mallinfo = mallinfo();
   
   /* Get total memory reserved by the system for malloc currently*/
   double reserved_mem = my_mallinfo.arena;
   
   /* Get all the memory currently allocated to user by malloc, etc. */
   double used_mem = my_mallinfo.hblkhd + my_mallinfo.usmblks +
      my_mallinfo.uordblks;
   
   /* Get memory not currently allocated to user but malloc controls */
   double free_mem = my_mallinfo.fsmblks + my_mallinfo.fordblks;
   
   /* Get number of items currently allocated */
   double number_allocated = my_mallinfo.ordblks + my_mallinfo.smblks;

   /* Record high-water mark for memory used. */
   s_max_memory = MathUtilities<double>::Max(s_max_memory,used_mem);   

   /* Print out concise malloc info line */
   os << used_mem/(1024.0*1024.0) << "MB ("
      << used_mem << ") in "
      << number_allocated << " allocs, "
      << reserved_mem/(1024.0*1024.0) << "MB reserved ("
      << free_mem << " unused)" << endl;
#endif
}

/*
*************************************************************************
*                                                                       *
* Records memory usage to user-defined event in TAU.  Note that if TAU  *
* is not included, this method does nothing.                            *
*                                                                       *
*************************************************************************
*/

void MemoryUtilities::recordMemoryInfo(double time) 
{

#ifdef HAVE_TAU
   /*
    * Access information from mallinfo
    */   
   struct mallinfo my_mallinfo = mallinfo();
   
   /* Get total memory reserved by the system for malloc currently*/
   double reserved_mem = my_mallinfo.arena;
   
   /* Get all the memory currently allocated to user by malloc, etc. */
   double used_mem = my_mallinfo.hblkhd + my_mallinfo.usmblks +
      my_mallinfo.uordblks;
   
   /* Get memory not currently allocated to user but malloc controls */
   double free_mem = my_mallinfo.fsmblks + my_mallinfo.fordblks;
   
   /* Get number of items currently allocated */
   double number_allocated = my_mallinfo.ordblks + my_mallinfo.smblks;


   /* These vars are unused now but we may use them in the future */
   NULL_USE(reserved_mem);
   NULL_USE(free_mem);
   NULL_USE(number_allocated);

   /* Record high-water mark for memory used. */
   s_max_memory = Utilities::dmax(s_max_memory,used_mem);   

   /*
    * Record "used_mem" in MB to tau event.
    */
   TAU_PROFILE_STMT(ue.TriggerEvent(used_mem/(1024.0*1024.0)));
   
#endif
}

/*
*************************************************************************
*                                                                       *
* Prints maximum memory used (i.e. high-water mark).  The max is        *
* determined each time the "printMemoryInfo" or "recordMemoryInfo"      *
* functions are called.                                                 *
*                                                                       *
*************************************************************************
*/
void MemoryUtilities::printMaxMemory(ostream& os) 
{
   /*
    * Step through all nodes (>0) and send max memory to processor 0,
    * which subsequently writes it out.
    */
   int maxmem = 0;
   int len = 1;
   for (int p = 0; p < MPI::getNodes(); p++) {
      if (MPI::getRank() == p) {
         maxmem = (int)s_max_memory;
         MPI::send(&maxmem, len, 0, false);
      }
      if (MPI::getRank() == 0) {
         MPI::recv(&maxmem, len, p, false);
      }
      os << "Maximum memory used on processor " << p
         << ": " << maxmem/(1024.*1024.) << " MB" << endl;
   }

}


}



