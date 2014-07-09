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
// File:	XDRStream.cc
// Package:	toolbox
// 
// 
// 
// Description:	Stream class that converts into XDR for portable communication
//

#include "toolbox/stream/XDRStream.h"

#include <string>
using namespace std;
#include "toolbox/base/Utilities.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "toolbox/stream/XDRStream.I"
#endif

namespace toolbox {



/*
*************************************************************************
*									*
* Some macros to simplify the implementation for the XDR streams	*
*									*
*************************************************************************
*/

#ifdef HAVE_XDR

#define XDR_PACK_OPAQUE(m_data,m_size)					\
   ((d_xdr_stream->x_op != XDR_ENCODE) ||				\
    (!xdr_opaque(d_xdr_stream, (caddr_t) m_data, m_size)))
#define XDR_UNPACK_OPAQUE(m_data,m_size)				\
   ((d_xdr_stream->x_op != XDR_DECODE) ||				\
    (!xdr_opaque(d_xdr_stream, (caddr_t) m_data, m_size)))
#define XDR_PACK_VECTOR(m_data,m_size,m_type)				\
   ((d_xdr_stream->x_op != XDR_ENCODE) ||				\
    (!xdr_vector(d_xdr_stream, (char *) m_data, m_size,			\
                 sizeof(m_type), (xdrproc_t) xdr_##m_type)))
#define XDR_UNPACK_VECTOR(m_data,m_size,m_type)				\
   ((d_xdr_stream->x_op != XDR_DECODE) ||				\
    (!xdr_vector(d_xdr_stream, (char *) m_data, m_size,			\
                 sizeof(m_type), (xdrproc_t) xdr_##m_type)))

#else

#define XDR_PACK_OPAQUE(m_data,m_size) 0
#define XDR_UNPACK_OPAQUE(m_data,m_size) 0
#define XDR_PACK_VECTOR(m_data,m_size,m_type) 0
#define XDR_UNPACK_VECTOR(m_data,m_size,m_type) 0

#endif


/*
*************************************************************************
*									*
* The virtual destructor for XDRStream does nothing.		*
*									*
*************************************************************************
*/

XDRStream::~XDRStream()
{
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for booleans.  Note that since	*
* the boolean representation is non-standard, boolean arrays are copied	*
* into character arrays and then packed using the character routines.	*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const bool& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(bool& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const bool *data, const int n)
{
   char *flags = new char[n];
   for (int i = 0; i < n; i++) {
      flags[i] = (data[i] ? 1 : 0);
   }
   if (XDR_PACK_OPAQUE(flags, n)) {
      TBOX_ERROR("XDRStream: Error in encoding bool...\n");
   }
   delete [] flags;
}

void XDRStream::unpack(bool *data, const int n)
{
   char *flags = new char[n];
   if (XDR_UNPACK_OPAQUE(flags, n)) {
      TBOX_ERROR("XDRStream: Error in decoding bool...\n");
   }
   for (int i = 0; i < n; i++) {
      data[i] = (flags[i] ? true : false);
   }
   delete [] flags;
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for characters			*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const char& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(char& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const char *data, const int n)
{
   if (XDR_PACK_OPAQUE(data, n)) {
      TBOX_ERROR("XDRStream: Error in encoding char...\n");
   }
}

void XDRStream::unpack(char *data, const int n)
{
   if (XDR_UNPACK_OPAQUE(data, n)) {
      TBOX_ERROR("XDRStream: Error in decoding char...\n");
   }
}

void XDRStream::writeString(const char *data)
{
#ifdef HAVE_XDR
   if (!xdr_string(d_xdr_stream, (char **) &data, strlen(data))) {
      TBOX_ERROR("XDRStream: Error in writing string...\n");
   }
#endif

}

#if 0
/*
*************************************************************************
*									*
* Packing and unpacking member functions for double complex		*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const dcomplex& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(dcomplex& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const dcomplex *data, const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
   if (XDR_PACK_VECTOR((double *) data, 2*n, double)) {
      TBOX_ERROR("XDRStream: Error in encoding double complex...\n");
   }
}

void XDRStream::unpack(dcomplex *data, const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
   if (XDR_UNPACK_VECTOR((double *) data, 2*n, double)) {
      TBOX_ERROR("XDRStream: Error in decoding double complex...\n");
   }
}
#endif

/*
*************************************************************************
*									*
* Packing and unpacking member functions for doubles			*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const double& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(double& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const double *data, const int n)
{
   if (XDR_PACK_VECTOR(data, n, double)) {
      TBOX_ERROR("XDRStream: Error in encoding double...\n");
   }
}

void XDRStream::unpack(double *data, const int n)
{
   if (XDR_UNPACK_VECTOR(data, n, double)) {
      TBOX_ERROR("XDRStream: Error in decoding double...\n");
   }
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for floats			*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const float& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(float& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const float *data, const int n)
{
   if (XDR_PACK_VECTOR(data, n, float)) {
      TBOX_ERROR("XDRStream: Error in encoding float...\n");
   }
}

void XDRStream::unpack(float *data, const int n)
{
   if (XDR_UNPACK_VECTOR(data, n, float)) {
      TBOX_ERROR("XDRStream: Error in decoding float...\n");
   }
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for integers			*
*									*
*************************************************************************
*/

AbstractStream& XDRStream::operator<<(const int& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& XDRStream::operator>>(int& data)
{
   unpack(&data, 1);
   return(*this);
}

void XDRStream::pack(const int *data, const int n)
{
   if (XDR_PACK_VECTOR(data, n, int)) {
      TBOX_ERROR("XDRStream: Error in encoding integer...\n");
   }
}

void XDRStream::unpack(int *data, const int n)
{
   if (XDR_UNPACK_VECTOR(data, n, int)) {
      TBOX_ERROR("XDRStream: Error in decoding integer...\n");
   }
}


}


