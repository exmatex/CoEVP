//
// File:	XDRStream.h
// Package:	toolbox
// 
// 
// 
// Description:	Stream class that converts into XDR for portable communication
//

#ifndef included_toolbox_XDRStream
#define included_toolbox_XDRStream


#ifdef HAVE_XDR

#ifndef included_sys_types
#define included_sys_types
extern "C" {
#include <sys/types.h>
}

#endif

#ifndef included_rpc_types
#define included_rpc_types
extern "C" {
#include <rpc/types.h>
}
#endif
#ifndef included_rpc_xdr
#define included_rpc_xdr
#ifndef __CYGWIN__
extern "C" {
#include <rpc/xdr.h>
}
#else
#define included_rpc_xdr
// Cygwin header file hack, xdr.h is not have ANSI C headers so 
// C++ thinks methods take no arguments
extern "C" {

#ifndef included_stdio
#define included_stdio
#include <stdio.h>
#endif

typedef	bool_t (*xdrproc_t)();
enum xdr_op {
	XDR_ENCODE=0,
	XDR_DECODE=1,
	XDR_FREE=2
};

typedef struct {
	enum xdr_op	x_op;		/* operation; fast additional param */
	struct xdr_ops {
		bool_t	(*x_getlong)();	/* get a long from underlying stream */
		bool_t	(*x_putlong)();	/* put a long to " */
		bool_t	(*x_getbytes)();/* get some bytes from " */
		bool_t	(*x_putbytes)();/* put some bytes to " */
		u_int	(*x_getpostn)();/* returns bytes off from beginning */
		bool_t  (*x_setpostn)();/* lets you reposition the stream */
		long *	(*x_inline)();	/* buf quick ptr to buffered data */
		void	(*x_destroy)();	/* free privates of this xdr_stream */
	} *x_ops;
	caddr_t 	x_public;	/* users' data */
	caddr_t		x_private;	/* pointer to private data */
	caddr_t 	x_base;		/* private used for position info */
	int		x_handy;	/* extra private word */
} XDR;

   extern bool_t xdr_opaque      (XDR *, caddr_t, u_int);
   extern bool_t xdr_string      (XDR *, char **, u_int);
   extern bool_t xdr_vector      (XDR *, char *, u_int, u_int, xdrproc_t);
   extern void   xdrmem_create   (XDR *, char *, u_int, enum xdr_op);
   extern void   xdrstdio_create (XDR *, FILE *, enum xdr_op);
   extern bool_t xdr_int	 (XDR *, int *);
   extern bool_t xdr_float 	 (XDR *, float *);
   extern bool_t xdr_double	 (XDR *, double *);
}
#endif
#endif

#else

typedef struct {
} XDR;

#endif 
// #ifdef HAVE_XDR

#ifndef included_toolbox_AbstractStream
#include "toolbox/stream/AbstractStream.h"
#endif

namespace toolbox {

/*!
 * @brief XDRStream provides pack/unpack operations that convert
 * between machine representation and XDR representation for portable
 * interprocessor communication.  Note that unless XDR translation is
 * required for portability across machine architectures, it would be
 * faster to copy data directly.
 *
 * The appropriate XDR stream must be set by setXDRStream() before any
 * of the packing or unpacking calls are made.  All of the packing and
 * unpacking operations have been defined here for the XDR operations.
 *
 * @see toolbox::AbstractStream
 */

class XDRStream : public AbstractStream
{
public:
   /*!
    * Standard default constructor for the XDR stream.
    */
   XDRStream();

   /*!
    * Set the XDR stream to/from which data will be written/read.  The
    * XDR stream must be set before calling any packing/unpacking operation.
    */
   void setXDRStream(XDR *xdrs);

   /*!
    * Virtual destructor for the XDR stream.
    */
   virtual ~XDRStream();

   /*!
    * @name Boolean Stream Primitives
    * Pack and unpack booleans into and out of the XDR data stream.
    */
   //@{
   /// Pack a single bool into the XDR data stream.
   virtual AbstractStream& operator<<(const bool& data);
   /// Remove a single bool from the XDR data stream.
   virtual AbstractStream& operator>>(bool& data);
   /// Pack an array of bools into the XDR data stream.
   virtual void pack(const bool* data, const int n = 1);
   /// Remove an array of bools from the XDR data stream.
   virtual void unpack(bool* data, const int n = 1);
   //@}

   /*!
    * @name Character Stream Primitives
    * Pack and unpack chars into and out of the XDR data stream.
    */
   //@{
   /// Pack a single char into the XDR data stream.
   virtual AbstractStream& operator<<(const char& data);
   /// Remove a single char from the XDR data stream.
   virtual AbstractStream& operator>>(char& data);
   /// Pack an array of chars into the XDR data stream.
   virtual void pack(const char* data, const int n = 1);
   /// Remove an array of chars from the XDR data stream.
   virtual void unpack(char* data, const int n = 1);
   /// Write a string into the XDR data stream.
   virtual void writeString(const char *data);
   //@}

#if 0
   /*!
    * @name Double Complex Stream Primitives
    * Pack and unpack double complex into and out of the message stream.
    */
   //@{
   /// Pack a single double complex into the message stream.
   virtual AbstractStream& operator<<(const dcomplex& data);
   /// Remove a single double complex from the message stream.
   virtual AbstractStream& operator>>(dcomplex& data);
   /// Pack an array of double complex into the message stream.
   virtual void pack(const dcomplex *data, const int n = 1);
   /// Remove an array of double complex from the message stream.
   virtual void unpack(dcomplex *data, const int n = 1);
   //@}
#endif

   /*!
    * @name Double Stream Primitives
    * Pack and unpack doubles into and out of the XDR data stream.
    */
   //@{
   /// Pack a single double into the XDR data stream.
   virtual AbstractStream& operator<<(const double& data);
   /// Remove a single double from the XDR data stream.
   virtual AbstractStream& operator>>(double& data);
   /// Pack an array of doubles into the XDR data stream.
   virtual void pack(const double* data, const int n = 1);
   /// Remove an array of doubles from the XDR data stream.
   virtual void unpack(double* data, const int n = 1);
   //@}

   /*!
    * @name Float Stream Primitives
    * Pack and unpack floats into and out of the XDR data stream.
    */
   //@{
   /// Pack a single float into the XDR data stream.
   virtual AbstractStream& operator<<(const float& data);
   /// Remove a single float from the XDR data stream.
   virtual AbstractStream& operator>>(float& data);
   /// Pack an array of floats into the XDR data stream.
   virtual void pack(const float* data, const int n = 1);
   /// Remove an array of floats from the XDR data stream.
   virtual void unpack(float* data, const int n = 1);
   //@}

   /*!
    * @name Integer Stream Primitives
    * Pack and unpack integers into and out of the XDR data stream.
    */
   //@{
   /// Pack a single integer into the XDR data stream.
   virtual AbstractStream& operator<<(const int& data);
   /// Remove a single integer from the XDR data stream.
   virtual AbstractStream& operator>>(int& data);
   /// Pack an array of integers into the XDR data stream.
   virtual void pack(const int* data, const int n = 1);
   /// Remove an array of integers from the XDR data stream.
   virtual void unpack(int* data, const int n = 1);
   //@}

private:
   XDR* d_xdr_stream;
};

}

#ifndef DEBUG_NO_INLINE
#include "toolbox/stream/XDRStream.I"
#endif

#endif
