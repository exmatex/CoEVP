// Copyright (c) 2002, A.H. van den Boogaard
//
// Test program for extended tensor functions

#include "xtensor.h"

using std::setw;
using std::cout;
using std::endl;
using std::setprecision;

int main()
{
   cout << setprecision(5);

   double A[36];
 
   // just some funny initialization
   for ( int i=0; i<6; i++ )
   {
      for ( int j=0; j<6; j++ )
         A[i+j*6] = (i+1)*(j+1) + (i+1)*i/double(j+1) + j/double(i+3);
 
      A[i+i*6] += 0.01;
   }

   Tensor4DSym T4D;
   Tensor4LSym T4L;

   T4D.putFortranMatrix( A ); // only lower left part used.
   T4L.putFortranMatrix( A );

   Tensor2Sym x, b;

   for ( int i=1; i<4; i++ )
      for ( int j=1; j<=i; j++ )
         b(i,j) = 1.0;

   x = solve( T4D, b );

   cout << "solve( Tensor4DSym, Tensor2Sym )" << endl;
   cout << "A" << endl << setw(12) << T4D << endl;
   cout << "b" << endl << setw(12) << b << endl;
   cout << "A:x=b -> x=" << endl << setw(12) << x << endl;
   cout << "residue A:x-b" << endl << setw(12) << (T4D*x-b) << endl;
 
   x = solve( T4L, b );

   cout << "solve( Tensor4LSym, Tensor2Sym )" << endl;
   cout << "A" << endl << setw(12) << T4L << endl;
   cout << "b" << endl << setw(12) << b << endl;
   cout << "A:x=b -> x=" << endl << setw(12) << x << endl;
   cout << "residue A:x-b" << endl << setw(12) << (T4L*x-b) << endl;

   Tensor4DSym T4Dinv( inv(T4D) );
   Tensor4LSym T4Linv( inv(T4L) );

   x = T4Dinv*b;
   cout << "T4Dinv*b" << endl << setw(12) << x << endl;

   x = T4Linv*b;
   cout << "T4Linv*b" << endl << setw(12) << x << endl;

   Tensor2Gen F(1);
   F(1,3) =  0.1;
   F(3,1) = -0.1;

   Tensor4DSym pfD;
   Tensor4LSym pfL;

   pfD = push_forward( T4D, F );
   pfL = push_forward( T4L, F );

   cout << "push forward pdD" << endl << setw(12) << pfD << endl;
   cout << "push forward pdL" << endl << setw(12) << pfL << endl;

   Tensor2Gen W(0), R;
   W(1,2) =  0.1;
   W(2,1) = -0.1;

   R = expW( W );

   cout << "Spin tensor" << endl << setw(12) << W << endl;
   cout << "Rotation tensor" << endl << setw(12) << R << endl;

//   getchar();  // for use with ms-windows

   return 0;
}
