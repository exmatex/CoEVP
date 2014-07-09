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
// File:	MessageStream.cc
// Package:	toolbox
// 
// 
// 
// Description:	Fixed-size message buffer used in interprocessor communication
//

#include "toolbox/stream/MessageStream.h"
#include "toolbox/base/Utilities.h"

#include <string.h>

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "toolbox/stream/MessageStream.I"
#endif

namespace toolbox {


bool MessageStream::s_use_xdr_translation = false;

/*
*************************************************************************
*									*
* The constructor and destructor for MessageStream.  The first	*
* form of the constructor uses the class-wide flag to determine whether	*
* XDR will be used.  The second constructor uses the use_xdr boolean	*
* argument.  If XDR translation is to be used, then create an XDR	*
* stream object	and bind it to the XDR manager object.			*
*									*
*************************************************************************
*/

MessageStream::MessageStream(const int bytes, const StreamMode mode)
{
   d_buffer_size  = bytes;
   d_current_size = 0;
   d_buffer_index = 0;
   d_use_xdr      = s_use_xdr_translation;
   d_buffer       = new char[d_buffer_size];

#ifdef HAVE_XDR
   if (d_use_xdr) {
      xdr_op xop = ((mode==MessageStream::Read) ? XDR_DECODE : XDR_ENCODE);
      xdrmem_create(&d_xdr_stream, (caddr_t) d_buffer, d_buffer_size, xop);
      d_xdr_manager.setXDRStream(&d_xdr_stream);
   }
#endif

}

MessageStream::MessageStream(const int bytes,
                                      const StreamMode mode,
                                      const bool use_xdr)
{
   d_buffer_size  = bytes;
   d_current_size = 0;
   d_buffer_index = 0;
   d_use_xdr      = use_xdr;
   d_buffer       = new char[d_buffer_size];

#ifdef HAVE_XDR
   if (d_use_xdr) {
      xdr_op xop = ((mode==MessageStream::Read) ? XDR_DECODE : XDR_ENCODE);
      xdrmem_create(&d_xdr_stream, (caddr_t) d_buffer, d_buffer_size, xop);
      d_xdr_manager.setXDRStream(&d_xdr_stream);
   }
#endif

}

MessageStream::~MessageStream()
{
#ifdef HAVE_XDR
   if (d_use_xdr) {
#ifndef LACKS_PROPER_XDR_HEADER
      xdr_destroy(&d_xdr_stream);
#else
      if (d_xdr_stream.x_ops->x_destroy)
         (*(void(*)(XDR*))(d_xdr_stream.x_ops->x_destroy))(&d_xdr_stream);
#endif
   }
#endif

   delete [] d_buffer;
}

/*
*************************************************************************
*									*
* Print out class data if an exception is thrown.			*
*									*
*************************************************************************
*/

void MessageStream::printClassData(ostream& os) const
{
   os << "Maximum buffer size = " << d_buffer_size << endl;
   os << "Current buffer size = " << d_current_size << endl;
   os << "Current buffer index = " << d_buffer_index << endl;
   os << "Pointer to buffer data = " << (void *) d_buffer << endl;
   os << "Using XDR translation = " << (d_use_xdr ? "true" : "false") << endl;
}

/*
*************************************************************************
*									*
* Packing/unpacking helper functions and macros.  The member function	*
* getPointerAndAdvanceCursor() returns a pointer to buffer space and	*
* advances internal pointers to reflect the allocated buffers space.	*
* The two macros given below simplify packing and unpacking for the	*
* numerous member functions below.					*
*									*
*************************************************************************
*/

void *MessageStream::getPointerAndAdvanceCursor(const int bytes)
{
   void *ptr = &d_buffer[d_buffer_index];
   d_buffer_index += bytes;
   if (d_buffer_index > d_current_size) {
      d_current_size = d_buffer_index;
      if (d_buffer_index > d_buffer_size) {
         TBOX_ERROR("MessageStream: Stream overrun of buffer...\n");
      }
   }
   return(ptr);
}

/*
*************************************************************************
*									*
* The following macros are used by all of the standard data types	*]
* except bool for packing and unpacking the data buffer.		*
*									*
*************************************************************************
*/

#ifdef HAVE_XDR

#define PACK(m_data,m_size,m_bytes) do {				\
   void *ptr = getPointerAndAdvanceCursor(m_bytes);			\
   if (d_use_xdr) {							\
      d_xdr_manager.pack(m_data, m_size);				\
   } else {								\
      memcpy(ptr, (void *) m_data, m_bytes);				\
   }                                                                    \
   } while (0)

#define UNPACK(m_data,m_size,m_bytes) do {				\
   void *ptr = getPointerAndAdvanceCursor(m_bytes);			\
   if (d_use_xdr) {							\
      d_xdr_manager.unpack(m_data, m_size);				\
   } else {								\
      memcpy((void *) m_data, ptr, m_bytes);				\
   }                                                                    \
   } while (0)

#else

#define PACK(m_data,m_size,m_bytes) do {				\
   void *ptr = getPointerAndAdvanceCursor(m_bytes);			\
   memcpy(ptr, (void *) m_data, m_bytes);                               \
   } while (0)
   
#define UNPACK(m_data,m_size,m_bytes) do {				\
   void *ptr = getPointerAndAdvanceCursor(m_bytes);			\
   memcpy((void *) m_data, ptr, m_bytes);                               \
   } while (0)

#endif


/*
*************************************************************************
*									*
* Packing and unpacking member functions for booleans.  Note that since	*
* the boolean representation is non-standard, boolean arrays are copied	*
* either using XDR or by converting into character arrays.		*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const bool& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(bool& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const bool *data, const int n)
{
   const int bytes = AbstractStream::sizeofBool(n);
   void *ptr = getPointerAndAdvanceCursor(bytes);
   if (d_use_xdr) {
#ifdef HAVE_XDR
      d_xdr_manager.pack(data, n);
#endif
   } else {
      char *c_ptr = (char *) ptr;
      for (int i = 0; i < n; i++) {
         c_ptr[i] = (data[i] ? 1 : 0);
      }
   }
}

void MessageStream::unpack(bool *data, const int n)
{
   const int bytes = AbstractStream::sizeofBool(n);
   void *ptr = getPointerAndAdvanceCursor(bytes);
   if (d_use_xdr) {
#ifdef HAVE_XDR
      d_xdr_manager.unpack(data, 1);
#endif
   } else {
      const char *c_ptr = (const char *) ptr;
      for (int i = 0; i < n; i++) {
         data[i] = (c_ptr[i] ? true : false);
      }
   }
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for characters			*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const char& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(char& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const char *data, const int n)
{
   const int bytes = AbstractStream::sizeofChar(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(char *data, const int n)
{
   const int bytes = AbstractStream::sizeofChar(n);
   UNPACK(data, n, bytes);
}

#if 0
/*
*************************************************************************
*									*
* Packing and unpacking member functions for double complex		*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const dcomplex& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(dcomplex& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const dcomplex *data, const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
   const int bytes = AbstractStream::sizeofDoubleComplex(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(dcomplex *data, const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
   const int bytes = AbstractStream::sizeofDoubleComplex(n);
   UNPACK(data, n, bytes);
}
#endif

/*
*************************************************************************
*									*
* Packing and unpacking member functions for doubles			*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const double& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(double& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const double *data, const int n)
{
   const int bytes = AbstractStream::sizeofDouble(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(double *data, const int n)
{
   const int bytes = AbstractStream::sizeofDouble(n);
   UNPACK(data, n, bytes);
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for floats			*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const float& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(float& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const float *data, const int n)
{
   const int bytes = AbstractStream::sizeofFloat(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(float *data, const int n)
{
   const int bytes = AbstractStream::sizeofFloat(n);
   UNPACK(data, n, bytes);
}

/*
*************************************************************************
*									*
* Packing and unpacking member functions for integers			*
*									*
*************************************************************************
*/

AbstractStream& MessageStream::operator<<(const int& data)
{
   pack(&data, 1);
   return(*this);
}

AbstractStream& MessageStream::operator>>(int& data)
{
   unpack(&data, 1);
   return(*this);
}

void MessageStream::pack(const int *data, const int n)
{
   const int bytes = AbstractStream::sizeofInt(n);
   PACK(data, n, bytes);
}

void MessageStream::unpack(int *data, const int n)
{
   const int bytes = AbstractStream::sizeofInt(n);
   UNPACK(data, n, bytes);
}


}


