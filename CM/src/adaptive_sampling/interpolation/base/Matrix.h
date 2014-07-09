//
// File:        Matrix.h
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

#if !defined(included_krigalg_Matrix)
#define included_krigalg_Matrix

#ifndef included_iosfwd
#define included_iosfwd
#include <iosfwd>
using namespace std;
#endif

#ifndef included_mtl_matrix
#define included_mtl_matrix
#include <mtl/matrix.h>
#endif

#ifndef included_krigalg_Vector
#include "Vector.h"
#endif

namespace krigalg {

  typedef mtl::matrix<double>::type Matrix;

  //
  // addition
  //

  Matrix & operator+=(Matrix       & matrix1,
		      const Matrix & matrix2);
  Matrix operator+(const Matrix & matrix1,
		   const Matrix & matrix2);

  //
  // subtraction
  //
  
  Matrix & operator-=(Matrix       & matrix1,
		      const Matrix & matrix2);
  Matrix operator-(const Matrix & matrix1,
		   const Matrix & matrix2);

  //
  // matrix multiplication
  //

  Matrix mult(const Matrix & matrix1,
	      const Matrix & matrix2);
  Matrix mult(const Matrix & matrix1,
	      const Matrix & matrix2,
	      bool           matrix1Transpose,
	      bool           matrix2Transpose);
	      
  //
  // multiplication of Vector by Matrix
  //

  Vector mult(const Matrix & matrix,
	      const Vector & vector);
  Vector mult(const Matrix & matrix,
	      const Vector & vector,
	      bool           matrixTranspose);

  //
  // generate transpose of Matrix
  //

  Matrix transpose(const Matrix & matrix);

  //
  // generate identity matrix
  //

  Matrix identity(int size);

  //
  // compute inverse
  //

  std::pair<Matrix, bool> inverse(const Matrix & matrix);

  //
  // output
  //
  
  std::ostream & operator<<(std::ostream & outputStream,
			    const Matrix & matrix);

}

#ifndef DEBUG_NO_INLINE
#include "Matrix.I"
#endif // DEBUG_NO_INLINE

#endif // included_krigalg_Matrix
