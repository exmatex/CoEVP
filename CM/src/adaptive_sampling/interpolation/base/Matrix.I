/* DO-NOT-DELETE revisionify.begin() */
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
/* DO-NOT-DELETE revisionify.end() */
//
// File:        Matrix.I
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing matrix.
//
// $Id$
//
// $Log$
//

#ifndef included_mtl_mtl
#define included_mtl_mtl
#include <mtl/mtl.h>
#endif

#ifndef included_mtl_scaled2D
#define included_mtl_scaled2D
#include <mtl/scaled2D.h>
#endif

#ifndef included_cassert
#define included_cassert
#include <cassert>
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_mtl_blas_interface_h
#define included_mtl_blas_interface_h
#include <mtl/blas_interface.h>
#endif 

#ifdef DEBUG_NO_INLINE
#define inline
#endif // DEBUG_NO_INLINE

namespace krigalg {
    
    //
    // addition
    //
    
    inline Matrix & 
    operator+=(Matrix       & matrix1,
	       const Matrix & matrix2)
    {
      
      mtl::add(matrix2,
	       matrix1);
      
      return matrix1;
      
    }
  
    
    inline Matrix 
    operator+(const Matrix & matrix1,
	      const Matrix & matrix2)
    {
      
      Matrix result(matrix1.nrows(),
		    matrix1.ncols());
      mtl::copy(matrix1, result);
      
      return result += matrix2;
      
    }
    
    //
    // subtraction
    //
    
    inline Matrix &
    operator-=(Matrix       & matrix1,
	       const Matrix & matrix2)
    {
      
      Matrix temp(matrix2.nrows(),
		  matrix2.ncols());
      
      mtl::copy(matrix2, temp);
      mtl::scale(temp, -1.0);
      
      mtl::add(temp,
	       matrix1);
      
      return matrix1;
      
    }

    inline Matrix
    operator-(const Matrix & matrix1,
	      const Matrix & matrix2)
    {
      
      Matrix result(matrix1.nrows(),
		    matrix1.ncols());
      mtl::copy(matrix1, result);
      
      return result -= matrix2;
      
    }
    
    //
    // matrix multiplication (A B)
    //
    
    inline Matrix
    mult(const Matrix & matrix1,
	 const Matrix & matrix2)
    {
      
      //
      // firewalls
      //
      
      assert(matrix1.ncols() == matrix2.nrows());
      
      //
      // instantiate Matrix of the proper dimensions
      //
      
      Matrix result(matrix1.nrows(),
		    matrix2.ncols());
      
      //
      // use mtl for multiplication
      //
      
#if 0
      
      mtl::mult(matrix1,
		matrix2,
		result);
      
#endif 
      
      //
      // since fortran stores arrays in the column order use 
      // C^T = B^T A^T to compute C=A B; the resulting array 
      // is stored properly and does not require transpose
      //
      //#if 0
      const char   transa = 'n';
      const char   transb = 'n';
      const int    m      = matrix2.ncols();
      const int    n      = matrix1.nrows();
      const int    k      = matrix2.nrows();
      const double alpha  = 1.0;
      const double * da   = matrix2.data();
      const int      lda  = matrix2.ncols();
      const double * db   = matrix1.data();
      const int      ldb  = matrix1.ncols();
      const double beta   = 0.0;
      double *     dc     = result.data();
      const int    ldc    = result.ncols();
      
      dgemm_(&transa, 
	     &transb,
	     &m, 
	     &n, 
	     &k,
	     &alpha, 
	     da,  
	     &lda,
	     db, 
	     &ldb,
	     &beta,
	     dc,
	     &ldc);
      //#endif     
      //
      //
      //

      return result;
      
    }
    
    //
    // matrix multiplication (A B || A^T B || A B^T || A^T B^T)
    //
    
    inline Matrix
    mult(const Matrix & matrix1,
	 const Matrix & matrix2,
	 bool           matrix1Transpose,
	 bool           matrix2Transpose)
    {

      //
      // no firewalls on the shape of the input matrices applied here:
      // blas performs that check
      //

      //
      // infer the shape of the result
      //

      int numberRows = matrix1Transpose ? matrix1.ncols() : matrix1.nrows();
      int numberCols = matrix2Transpose ? matrix2.nrows() : matrix2.ncols();

      //
      // instantiate Matrix of the proper dimensions
      //
      
      Matrix result(numberRows,
		    numberCols);
      
      //
      // since fortran stores arrays in the column order use 
      // C^T = B^T A^T to compute C=A B; the resulting array 
      // is stored properly and does not require transpose
      //

      const char   transa = matrix2Transpose ? 't' : 'n';
      const char   transb = matrix1Transpose ? 't' : 'n';
      const int    m      = matrix2Transpose ? matrix2.nrows():matrix2.ncols();
      const int    n      = matrix1Transpose ? matrix1.ncols():matrix1.nrows();
      const int    k      = matrix2Transpose ? matrix2.ncols():matrix2.nrows();
      const double alpha  = 1.0;
      const double * da   = matrix2.data();
      const int      lda  = matrix2Transpose ? k : m;
      const double * db   = matrix1.data();
      const int      ldb  = matrix1Transpose ? n : k;
      const double beta   = 0.0;
      double *     dc     = result.data();
      const int    ldc    = result.ncols();
      
      dgemm_(&transa, 
	     &transb,
	     &m, 
	     &n, 
	     &k,
	     &alpha, 
	     da,  
	     &lda,
	     db, 
	     &ldb,
	     &beta,
	     dc,
	     &ldc);

      //
      //
      //

      return result;

    }

    //
    // multiplication of Vector by Matrix
    //
    
    inline Vector 
    mult(const Matrix & matrix,
	 const Vector & vector)
    {
      
      //
      // firewalls
      //
      
      assert(matrix.ncols() == vector.size());
      
      //
      // instantiate Vector of size matrix.major()
      //
      
      Vector result(matrix.nrows());
      
      //
      // use mtl for multiplication
      //
      
#if 0
      
      mtl::mult(matrix,
		vector,
		result);
      
#endif 
      
      //#if 0
      
      char trans = 't';
      int m    = matrix.ncols();
      int n    = matrix.nrows();
      double alpha      = 1.0;
      Matrix::value_type * A = 
	const_cast<Matrix::value_type *>(matrix.data()); 
      int lda           = matrix.ncols();
      int incx          = 1;
      double * X        = const_cast<double *>(&(vector[0]));
      double   beta     = 0.0;
      double * Y        = &(result[0]);
      
      dgemv_(&trans,
	     &m,
	     &n,
	     &alpha,
	     A,
	     &lda,
	     X,
	     &incx,
	     &beta,
	     Y,
	     &incx);
      
      //#endif 

      //
      //
      //

    return result;

    }

    inline Vector
    mult(const Matrix & matrix,
	 const Vector & vector,
	 bool           matrixTranspose)
    {

      //
      // determine the size of the resulting vector
      //

      const int resultSize = matrixTranspose ? matrix.ncols() : matrix.nrows();

      //
      // firewalls
      //
      
      if (matrixTranspose == false)
	assert(matrix.ncols() == vector.size());
      else
	assert(matrix.nrows() == vector.size());

      //
      // instantiate Vector of size matrix.major()
      //
      
      Vector result(resultSize);
      
      //
      // use blas for multiplication
      //
      
      char trans        = matrixTranspose ? 'n' : 't';
      int m             = matrix.ncols();
      int n             = matrix.nrows();
      double alpha      = 1.0;
      Matrix::value_type * A = 
	const_cast<Matrix::value_type *>(matrix.data()); 
      int lda           = matrix.ncols();
      int incx          = 1;
      double * X        = const_cast<double *>(&(vector[0]));
      double   beta     = 0.0;
      double * Y        = &(result[0]);
      
      dgemv_(&trans,
	     &m,
	     &n,
	     &alpha,
	     A,
	     &lda,
	     X,
	     &incx,
	     &beta,
	     Y,
	     &incx);
      

      //
      //
      //

      return result;

    }

    //
    // generate transpose of Matrix
    //
    
    inline Matrix 
    transpose(const Matrix & matrix)
    {
      
      Matrix result(matrix.ncols(),
		    matrix.nrows());
      mtl::transpose(matrix, 
		     result);
      
      return result;
      
    }

    //
    // generate identity matrix
    //
    
    inline Matrix 
    identity(int size)
    {
      
      Matrix result(size,
		    size);
      mtl::set_diagonal(result,
			1.0);
      return result;
      
    }

    //
    // output operator
    //
    
    inline std::ostream & 
    operator<<(std::ostream & outputStream,
	       const Matrix & matrix)
    {

      typedef mtl::matrix_traits<Matrix>::size_type Int;
      Int i,j;
      outputStream << matrix.nrows() << "x" << matrix.ncols() << std::endl;
      for (i=0; i < matrix.nrows(); ++i) {
	for (j=0; j < matrix.ncols(); ++j) {
	  outputStream << matrix(i,j);
	  if (j < matrix.ncols() - 1)
	    outputStream << "\t";
	}
	if (i < matrix.nrows() - 1)
	  outputStream << "\t" << std::endl;
	else
	  outputStream << std::endl;
      }

      return outputStream;

    }

}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif // DEBUG_NO_INLINE


