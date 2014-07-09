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
// File:  MPI.cc
// Package:  toolbox
// 
// 
// 
// Description:  Simple utility class for interfacing with MPI
//

#include "toolbox/parallel/MPI.h"

#include <stdlib.h>
#include <string.h>

#include "toolbox/base/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "toolbox/parallel/MPI.I"
#endif

#ifdef __INSURE__
/*
 * These are defined in mpich mpi.h and break the insure compile.
 * This may impact Globus in some way, at least from the comments
 * in the mpi.h header file.  Why mpich externs something that is 
 * not defined in the mpich is confusing and probably just broken.
 */
int MPICHX_TOPOLOGY_DEPTHS;
int MPICHX_TOPOLOGY_COLORS;
#endif



namespace toolbox {

MPI::comm MPI::s_communicator      = (MPI::comm) 0;
int      MPI::s_outgoing_messages = 0;
int      MPI::s_outgoing_bytes    = 0;
int      MPI::s_incoming_messages = 0;
int      MPI::s_incoming_bytes    = 0;
int      MPI::s_initialized       = 0; 

#if HAVE_MPI
MPI::comm MPI::commWorld = MPI_COMM_WORLD;
MPI::comm MPI::commNull = MPI_COMM_NULL;
#else
MPI::comm MPI::commWorld = 0;
MPI::comm MPI::commNull = -1;
#endif

bool MPI::s_call_abort_in_serial_instead_of_exit = true;

/*
**************************************************************************
*                                                                        *
* Abort the program.                                                     *
*                                                                        *
**************************************************************************
*/

void MPI::setCallAbortInSerialInsteadOfExit(bool flag)
{
   s_call_abort_in_serial_instead_of_exit = flag;
}

void MPI::abort()
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      MPI_Abort(s_communicator, -1);
   } else {
      if (s_call_abort_in_serial_instead_of_exit) {
         ::abort();
      } else {
         exit(-1);
      }
   }
#else
   if (s_call_abort_in_serial_instead_of_exit) {
      ::abort();
   } else {
      exit(-1);
   }
#endif

}

/*
**************************************************************************
*                                                                        *
* Initialize the static data in the MPI utility class.  This must be     *
* called after MPI_Init to ensure that the MPI_COMM_WORLD structure      *
* has been initialized.                                                  *
*                                                                        *
**************************************************************************
*/

void MPI::initialize()
{
   if (!s_initialized) {
      s_communicator      = MPI::commWorld;
      s_outgoing_messages = 0;
      s_outgoing_bytes    = 0;
      s_incoming_messages = 0;
      s_incoming_bytes    = 0;
      s_initialized = 1;
   }
}

/*
**************************************************************************
*                                                                        *
* Tree depth calculation for tracking the * number of message sends      *
* and receives and the number of bytes.                                  *
*                                                                        *
**************************************************************************
*/

int MPI::getTreeDepth()
{
   int depth = 0;
   const int nnodes = getNodes();
   while ((1 << depth) < nnodes) {
      depth++;
   }
   return(depth);
}

/*
**************************************************************************
*                                                                        *
* Perform a global barrier across all processors.                        *
*                                                                        *
**************************************************************************
*/

void MPI::barrier()
{
#ifdef HAVE_MPI
   (void) MPI_Barrier(s_communicator);
   const int tree = getTreeDepth();
   updateOutgoingStatistics(tree, 0);
   updateIncomingStatistics(tree, 0);
#endif
}
 
/*
**************************************************************************
*                                                                        *
* Sum reduction for a scalar double.                                     *
*                                                                        *
**************************************************************************
*/

double MPI::sumReduction(const double x)
{
   double recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      double send = x;
      MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE, MPI_SUM, s_communicator);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(double));
      updateIncomingStatistics(tree, sizeof(double));
   }
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for an array of doubles.                                 *
*                                                                        *
**************************************************************************
*/

void MPI::sumReduction(double *x, const int n)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      double *send = new double[n];
      memcpy(send, x, n*sizeof(double));
      MPI_Allreduce(send, x, n, MPI_DOUBLE, MPI_SUM, s_communicator);
      delete [] send;
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(double));
      updateIncomingStatistics(tree, n*sizeof(double));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for a scalar float.                                      *
*                                                                        *
**************************************************************************
*/

float MPI::sumReduction(const float x)
{
   float recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      float send = x;
      MPI_Allreduce(&send, &recv, 1, MPI_FLOAT, MPI_SUM, s_communicator);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(float));
      updateIncomingStatistics(tree, sizeof(float));
   }
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for an array of floats.                                  *
*                                                                        *
**************************************************************************
*/

void MPI::sumReduction(float *x, const int n)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      float *send = new float[n];
      memcpy(send, x, n*sizeof(float));
      MPI_Allreduce(send, x, n, MPI_FLOAT, MPI_SUM, s_communicator);
      delete [] send;
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(float));
      updateIncomingStatistics(tree, n*sizeof(float));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

#if 0
/*
**************************************************************************
*                                                                        *
* Sum reduction for a scalar dcomplex.                                   *
*                                                                        *
**************************************************************************
*/

dcomplex MPI::sumReduction(const dcomplex x)
{
   dcomplex recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      double xreal[2];
      double send[2];
      send[0] = xreal[0] = real(x); 
      send[1] = xreal[1] = imag(x);
      MPI_Allreduce(send, xreal, 2, MPI_DOUBLE, MPI_SUM, s_communicator);
      recv = dcomplex(xreal[0], xreal[1]);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, 2*sizeof(double));
      updateIncomingStatistics(tree, 2*sizeof(double));
   }
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for an array of dcomplex.                                *
*                                                                        *
**************************************************************************
*/

void MPI::sumReduction(dcomplex *x, const int n)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      int nrvals = 2*n;
      double *xreal = new double[nrvals];
      double *send = new double[nrvals];
      for (int i = 0; i < n; i++) {
         xreal[2*i]   = real(x[i]); 
         xreal[2*i+1] = imag(x[i]); 
      }
      memcpy(send, xreal, nrvals*sizeof(double));
      MPI_Allreduce(send, xreal, nrvals, MPI_DOUBLE, MPI_SUM, s_communicator);
      for (int j = 0; j < n; j++) {
         x[j] = dcomplex(xreal[2*j], xreal[2*j+1]);
      }
      delete [] send;
      delete [] xreal;
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(double));
      updateIncomingStatistics(tree, n*sizeof(double));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}
#endif

/*
**************************************************************************
*                                                                        *
* Sum reduction for a scalar integer.                                    *
*                                                                        *
**************************************************************************
*/

int MPI::sumReduction(const int x)
{
   int recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      int send = x;
      MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_SUM, s_communicator);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(int));
      updateIncomingStatistics(tree, sizeof(int));
   }
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Sum reduction for an array of integers.                                z
*                                                                        *
**************************************************************************
*/

void MPI::sumReduction(int *x, const int n)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      int *send = new int[n];
      memcpy(send, x, n*sizeof(int));
      MPI_Allreduce(send, x, n, MPI_INT, MPI_SUM, s_communicator);
      delete [] send;
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(int));
      updateIncomingStatistics(tree, n*sizeof(int));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Min reduction for a scalar double.                                     *
*                                                                        *
**************************************************************************
*/

double MPI::minReduction(const double x, int *rank_of_min)
{
   double rval = x;

   /*
    * If a rank_of_min argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_min != NULL) {
      *rank_of_min = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         double send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_DOUBLE, MPI_MIN, s_communicator);
      } else {
         DoubleIntStruct recv;
         DoubleIntStruct send;
         send.d = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_DOUBLE_INT,
                       MPI_MINLOC,
                       s_communicator);
         rval = recv.d;
         *rank_of_min = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(double));
      updateIncomingStatistics(tree, sizeof(double));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Min reduction for an array of doubles.                                 *
*                                                                        *
**************************************************************************
*/

void MPI::minReduction(double *x, const int n, int *rank_of_min)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         double *send = new double[n];
         memcpy(send, x, n*sizeof(double));
         MPI_Allreduce(send, x, n, MPI_DOUBLE, MPI_MIN, s_communicator);
         delete [] send;
      }
      else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_DOUBLE_INT,
                       MPI_MINLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].d;
            rank_of_min[i] = send[i].i;
         }
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(double));
      updateIncomingStatistics(tree, n*sizeof(double));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Min reduction for a scalar float.                                      *
*                                                                        *
**************************************************************************
*/

float MPI::minReduction(const float x, int *rank_of_min)
{
   float rval = x;

   /*
    * If a rank_of_min argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_min != NULL) {
      *rank_of_min = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         float send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_FLOAT, MPI_MIN, s_communicator);
      }
      else {
         FloatIntStruct recv;
         FloatIntStruct send;
         send.f = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_FLOAT_INT,
                       MPI_MINLOC,
                       s_communicator);
         rval = recv.f;
         *rank_of_min = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(float));
      updateIncomingStatistics(tree, sizeof(float));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Min reduction for an array of floats.                                  *
*                                                                        *
**************************************************************************
*/

void MPI::minReduction(float *x, const int n, int *rank_of_min)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         float *send = new float[n];
         memcpy(send, x, n*sizeof(float));
         MPI_Allreduce(send, x, n, MPI_FLOAT, MPI_MIN, s_communicator);
         delete [] send;
      }
      else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_FLOAT_INT,
                       MPI_MINLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].f;
            rank_of_min[i] = send[i].i;
         }
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(float));
      updateIncomingStatistics(tree, n*sizeof(float));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Min reduction for a scalar integer.                                    *
*                                                                        *
**************************************************************************
*/

int MPI::minReduction(const int x, int *rank_of_min)
{
   int rval = x;

   /*
    * If a rank_of_min argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_min != NULL) {
      *rank_of_min = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         int send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_INT, MPI_MIN, s_communicator);
      }
      else {
         IntIntStruct recv;
         IntIntStruct send;
         send.j = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_2INT,
                       MPI_MINLOC,
                       s_communicator);
         rval = recv.j;
         *rank_of_min = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(int));
      updateIncomingStatistics(tree, sizeof(int));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Min reduction for an array of integers.                                *
*                                                                        *
**************************************************************************
*/

void MPI::minReduction(int *x, const int n, int *rank_of_min)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_min == NULL ) {
         int *send = new int[n];
         memcpy(send, x, n*sizeof(int));
         MPI_Allreduce(send, x, n, MPI_INT, MPI_MIN, s_communicator);
         delete [] send;
      }
      else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_2INT,
                       MPI_MINLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            rank_of_min[i] = send[i].i;
         }
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(int));
      updateIncomingStatistics(tree, n*sizeof(int));
   }
#else
   NULL_USE(x);
   NULL_USE(n);

#endif
}

/*
**************************************************************************
*                                                                        *
* Max reduction for a scalar double.                                     *
*                                                                        *
**************************************************************************
*/

double MPI::maxReduction(const double x, int *rank_of_max)
{
   double rval = x;

   /*
    * If a rank_of_max argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_max != NULL) {
      *rank_of_max = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         double send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_DOUBLE, MPI_MAX, s_communicator);
      } else {
         DoubleIntStruct recv;
         DoubleIntStruct send;
         send.d = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_DOUBLE_INT,
                       MPI_MAXLOC,
                       s_communicator);
         rval = recv.d;
         *rank_of_max = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(double));
      updateIncomingStatistics(tree, sizeof(double));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Max reduction for an array of doubles.                                 *
*                                                                        *
**************************************************************************
*/

void MPI::maxReduction(double *x, const int n, int *rank_of_max)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         double *send = new double[n];
         memcpy(send, x, n*sizeof(double));
         MPI_Allreduce(send, x, n, MPI_DOUBLE, MPI_MAX, s_communicator);
         delete [] send;
      }
      else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_DOUBLE_INT,
                       MPI_MAXLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].d;
            rank_of_max[i] = send[i].i;
         }
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(double));
      updateIncomingStatistics(tree, n*sizeof(double));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Max reduction for a scalar float.                                      *
*                                                                        *
**************************************************************************
*/

float MPI::maxReduction(const float x, int *rank_of_max)
{
   float rval = x;

   /*
    * If a rank_of_max argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_max != NULL) {
      *rank_of_max = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         float send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_FLOAT, MPI_MAX, s_communicator);
      }
      else {
         FloatIntStruct recv;
         FloatIntStruct send;
         send.f = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_FLOAT_INT,
                       MPI_MAXLOC,
                       s_communicator);
         rval = recv.f;
         *rank_of_max = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(float));
      updateIncomingStatistics(tree, sizeof(float));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Max reduction for an array of floats.                                  *
*                                                                        *
**************************************************************************
*/

void MPI::maxReduction(float *x, const int n, int *rank_of_max)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         float *send = new float[n];
         memcpy(send, x, n*sizeof(float));
         MPI_Allreduce(send, x, n, MPI_FLOAT, MPI_MAX, s_communicator);
         delete [] send;
      }
      else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_FLOAT_INT,
                       MPI_MAXLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].f;
            rank_of_max[i] = send[i].i;
         }
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(float));
      updateIncomingStatistics(tree, n*sizeof(float));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* Max reduction for a scalar integer.                                    *
*                                                                        *
**************************************************************************
*/

int MPI::maxReduction(const int x, int *rank_of_max)
{
   int rval = x;

   /*
    * If a rank_of_max argument is provided, set it to the current
    * rank of the process.
    */
   if (rank_of_max != NULL) {
      *rank_of_max = getRank();
   }

#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         int send = x;
         MPI_Allreduce(&send, &rval, 1, MPI_INT, MPI_MAX, s_communicator);
      }
      else {
         IntIntStruct recv;
         IntIntStruct send;
         send.j = x;
         send.i = getRank();
         MPI_Allreduce(&send,
                       &recv,
                       1,
                       MPI_2INT,
                       MPI_MAXLOC,
                       s_communicator);
         rval = recv.j;
         *rank_of_max = recv.i;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(int));
      updateIncomingStatistics(tree, sizeof(int));
   }
#endif
   return(rval);
}

/*
**************************************************************************
*                                                                        *
* Max reduction for an array of integers.                                *
*                                                                        *
**************************************************************************
*/

void MPI::maxReduction(int *x, const int n, int *rank_of_max)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      if ( rank_of_max == NULL ) {
         int *send = new int[n];
         memcpy(send, x, n*sizeof(int));
         MPI_Allreduce(send, x, n, MPI_INT, MPI_MAX, s_communicator);
         delete [] send;
      }
      else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = getRank();
         }
         MPI_Allreduce(send,
                       recv,
                       n,
                       MPI_2INT,
                       MPI_MAXLOC,
                       s_communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            rank_of_max[i] = send[i].i;
         }
	 delete recv;
	 delete send;
      }
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, n*sizeof(int));
      updateIncomingStatistics(tree, n*sizeof(int));
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}

/*
**************************************************************************
*                                                                        *
* All-to-one sum reduction on integer array.                             *
*                                                                        *
**************************************************************************
*/

void MPI::allToOneSumReduction(int *x, const int n, const int root)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      int *send = new int[n];
      memcpy(send, x, n*sizeof(int));
      MPI_Reduce(send, x, n, MPI_INT, MPI_SUM, root, s_communicator);
      delete [] send;
      const int tree = getTreeDepth();

      /*
       * note: following probably isn't correct, since some processors
       * both send and receive, and some only send; this is of course
       * dependent on the MPI implementation and/or the hardware.
       */
      if (getRank() == root) {
         updateOutgoingStatistics(tree, n*sizeof(int));
      } else {
         updateIncomingStatistics(tree, n*sizeof(int));
      }
   }
#else
   NULL_USE(x);
   NULL_USE(n);
#endif
}


/*
**************************************************************************
*                                                                        *
* Broadcast for scalar integer.                                          *
*                                                                        *
**************************************************************************
*/

int MPI::bcast(const int x, const int root)
{
   int recv = x;
#ifdef HAVE_MPI
   if (getNodes() > 1) {
      (void) MPI_Bcast(&recv, 1, MPI_INT, root, s_communicator);
      const int tree = getTreeDepth();
      if (getRank() == root) {
         updateOutgoingStatistics(tree, sizeof(int));
      } else {
         updateIncomingStatistics(tree, sizeof(int));
      }
   }
#else
   NULL_USE(x);
   NULL_USE(root);
#endif
   return(recv);
}

/*
**************************************************************************
*                                                                        *
* Broadcast for integer array from root processor to all other procs.    *
*                                                                        *
**************************************************************************
*/

void MPI::bcast(int *x, int &length, const int root)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {

      (void) MPI_Bcast((void*)x, length, MPI_INT, root, s_communicator);
      const int tree = getTreeDepth();
      if (getRank() == root) {
         updateOutgoingStatistics(tree, length*sizeof(int));
      } else {
         updateIncomingStatistics(tree, length*sizeof(int));
      }
   }
#else
   NULL_USE(x);
   NULL_USE(root);
#endif
}

/*
**************************************************************************
*                                                                        *
* Broadcast for char array from root processor to all other processors.  *
*                                                                        *
**************************************************************************
*/

void MPI::bcast(char *x, int &length, const int root)
{
#ifdef HAVE_MPI
   if (getNodes() > 1) {

      (void) MPI_Bcast((void*)x, length, MPI_BYTE, root, s_communicator);
      const int tree = getTreeDepth();
      if (getRank() == root) {
         updateOutgoingStatistics(tree, length*sizeof(int));
      } else {
         updateIncomingStatistics(tree, length*sizeof(int));
      }
   }
#else
   NULL_USE(x);
   NULL_USE(root);
#endif
}

/*
**************************************************************************
*                                                                        *
* Send integer array to another processor.                               *
*                                                                        *
**************************************************************************
*/

void MPI::send(const int *buf, 
                    const int length, 
                    const int receiving_proc_number, 
                    const bool send_length, 
                    int tag)
{
#ifdef HAVE_MPI
   tag = (tag >= 0) ? tag : 0;
   int size = length;
   if (send_length) {
      MPI_Send(&size, 1, MPI_INT, receiving_proc_number, tag, s_communicator);
      const int tree = getTreeDepth();
      updateOutgoingStatistics(tree, sizeof(int));
   }
   MPI_Send((void*)buf, 
            length, 
            MPI_INT, 
            receiving_proc_number, 
            tag, 
            s_communicator);
   const int tree = getTreeDepth();
   updateOutgoingStatistics(tree, length * sizeof(int));
#endif
}

/*
**************************************************************************
*                                                                        *
* Receive integer array from another processor.                          *
*                                                                        *
**************************************************************************
*/

void MPI::recv(int *buf, 
                    int &length, 
                    const int sending_proc_number,
		    const bool get_length, 
                    int tag)
{
#ifdef HAVE_MPI
   MPI_Status status;
   tag = (tag >= 0) ? tag : 0;
   if (get_length) {
      MPI_Recv(&length, 
               1, 
               MPI_INT, 
               sending_proc_number, 
               tag, 
               s_communicator,
               &status);
      const int tree = getTreeDepth();
      updateIncomingStatistics(tree, sizeof(int));
   }
   MPI_Recv((void*)buf, 
            length, 
            MPI_INT, 
            sending_proc_number, 
            tag, 
            s_communicator, 
            &status);
   const int tree = getTreeDepth();
   updateIncomingStatistics(tree, length*sizeof(int));
#endif
}


/*
*************************************************************************
*									*
* Send an array of number_bytes bytes from this processer to            *
* receiving_proc.                                                       *
* This call must be paired with a matching call to MPI::recvBytes. *
*									*
*************************************************************************
*/

void MPI::sendBytes(const void *buf, 
                         const int number_bytes, 
                         const int receiving_proc_number)
{
#ifdef HAVE_MPI

   MPI_Send((void*)buf, 
            number_bytes, 
            MPI_BYTE, 
            receiving_proc_number, 
            0, 
            s_communicator);
   const int tree = getTreeDepth();
   updateOutgoingStatistics(tree, number_bytes * sizeof(char));
#endif
}



/*
*************************************************************************
*									*
* Receive an array of bytes of max size number_bytes bytes from any     *
* processer.                                                            *
* This call must be paired with a matching call to MPI::sendBytes. *
*                                                                       *
* Returns the processor number of the sender.                           *
*									*
*************************************************************************
*/

int MPI::recvBytes(void *buf, 
                        int number_bytes) 
{
   int rval = 0;
#ifdef HAVE_MPI
   MPI_Status status;
   MPI_Recv(buf, 
            number_bytes, 
            MPI_BYTE, 
            MPI_ANY_SOURCE, 
            MPI_ANY_TAG, 
            s_communicator, 
            &status);

   const int tree = getTreeDepth();
   updateIncomingStatistics(tree, number_bytes * sizeof(char));
   rval = status.MPI_SOURCE;
#endif

   return rval;
}

/*
**************************************************************************
*                                                                        *
* All-to-all exchange of arrays of integers; each processor's            *
* array can be of a different length                                     *
*                                                                        *
**************************************************************************
*/

void MPI::allGather(
   const int *x_in, int size_in, int *x_out, int size_out)
{
#ifdef HAVE_MPI
   int* rcounts = (int*)NULL;
   int* disps = (int*)NULL;
   allGatherSetup(size_in, size_out, rcounts, disps);

   MPI_Allgatherv((void*)x_in, size_in, MPI_INT,
                  x_out, rcounts, disps, MPI_INT, s_communicator);

   if (rcounts) {
      delete [] rcounts;
   }
   if (disps) {
      delete [] disps;
   }
#else
   NULL_USE(x_in);
   NULL_USE(size_in);
   NULL_USE(x_out);
   NULL_USE(size_out);
#endif
}

/*
**************************************************************************
*                                                                        *
* all-to-all exchange of arrays of doubles; each processor's             *
* array can be of a different length                                     *
*                                                                        *
**************************************************************************
*/

void MPI::allGather(
   const double *x_in, int size_in, double *x_out, int size_out)
{
#ifdef HAVE_MPI
   int *rcounts = (int*)NULL;
   int *disps = (int*)NULL;
   allGatherSetup(size_in, size_out, rcounts, disps);

   MPI_Allgatherv((void*)x_in, size_in, MPI_DOUBLE,
                  x_out, rcounts, disps, MPI_DOUBLE, s_communicator);

   if (rcounts) {
      delete [] rcounts;
   }
   if (disps) {
      delete [] disps;
   }
#else
   NULL_USE(x_in);
   NULL_USE(size_in);
   NULL_USE(x_out);
   NULL_USE(size_out);
#endif
}

/*
*************************************************************************
*                                                                       *
* common setup funtion for all-to-all functions                         *
*                                                                       *
*************************************************************************
*/

void MPI::allGatherSetup(
   int size_in, int size_out, int *&rcounts, int *&disps)
{
#ifdef HAVE_MPI
   int np = MPI::getNodes();
   rcounts = new int[np];
   disps = new int[np];

   /* figure out where where each processor's input will be placed */
   MPI::allGather(size_in, rcounts);

   disps[0] = 0;
   for (int p = 1; p < np; ++p) {
      disps[p] = disps[p-1] + rcounts[p-1];
   }

   /* verify that the x_out array is the appropriate size! */
   int c = 0;
   for (int x = 0; x < np; ++x) {
      c += rcounts[x];
   }
   if (c != size_out) {
      TBOX_ERROR("MPI::allGatherSetup error..." 
                 << "\n   size_out =" << size_out << "appears to be incorrect; "
                 << "should be: " << c << endl);
   }
#else
   NULL_USE(size_in);
   NULL_USE(size_out);
#endif
}

/*
**************************************************************************
*                                                                        *
* All-to-all exchange of a double.                                       *
*                                                                        *
**************************************************************************
*/

void MPI::allGather(double x_in, double *x_out)
{
#ifdef HAVE_MPI
   MPI_Allgather(&x_in, 1, MPI_DOUBLE, x_out, 1, MPI_DOUBLE, s_communicator);
#else
   NULL_USE(x_in);
   NULL_USE(x_out);
#endif
}

/*
**************************************************************************
*                                                                        *
* All-to-all exchange of an integer.                                     *
*                                                                        *
* ************************************************************************
*/

void MPI::allGather(int x_in, int *x_out)
{
#ifdef HAVE_MPI
   MPI_Allgather(&x_in, 1, MPI_INT, x_out, 1, MPI_INT, s_communicator);
#else
   NULL_USE(x_in);
   NULL_USE(x_out);
#endif
}

}


