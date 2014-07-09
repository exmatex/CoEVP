// Copyright (c) 2001, A.H. van den Boogaard
//
// Test program for tensor class

#include "tensor.h"
#include <cmath>

using std::setw;
using std::cout;
using std::endl;
using std::setprecision;

int main()
{
   cout << setprecision(5);

   Tensor1    a, b;
   Tensor2Gen A, B, C;

   a(1) = 1; a(2) = 2; a(3) = 3;
   b(1) = 5; b(2) = 6; b(3) = 7;

   A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
   A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
   A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

   cout << "a" << endl << setw(12) << a << endl;
   cout << "b" << endl << setw(12) << b << endl;
   cout << "A" << endl << setw(12) << A << endl;

   B = dyadic(a,b);
   C = dyadic(b,a);

   cout << "B = dyadic(a,b)" << endl << setw(12) << B << endl;
   cout << "C = dyadic(b,a)" << endl << setw(12) << C << endl;

   cout << "sym(C) " << endl << setw(12) << sym(C) << endl;
   cout << "skew(C) " << endl << setw(12) << skew(C) << endl;
   cout << "sym(C) + skew(C) " << endl << setw(12) << (sym(C)+skew(C)) << endl;

   A = B*C;

   cout << "B*C" << endl << setw(12) << A << endl;
   cout << "trace(A) " << setw(12) <<  trace(A) << endl;

   cout << "dev(A) " << endl << setw(12) << dev(A) << endl;


   // Examples from computational mechanics

   Tensor2Gen F, Finv, R;
   Tensor2Sym U;

   F(1,1) =  2  ; F(1,2) =  0  ; F(1,3) = 0.5;
   F(2,1) = -0.3; F(2,2) =  1  ; F(2,3) = 0.5;
   F(3,1) =  0  ; F(3,2) = -0.2; F(3,3) = 0.8;

   Finv = inv(F);

   cout << "F" << endl << setw(12) << F << endl;
   cout << "inv(F)" << endl << setw(12) << Finv << endl;
   cout << "F*inv(F)" << endl << setw(12) << F*Finv << endl;

   polarRight( F, &R, &U );

   cout << "polar decomposition F=RU" << endl;
   cout << "R" << endl << setw(12) << R << endl;
   cout << "U" << endl << setw(12) << U << endl;
   cout << "R.U" << endl << setw(12) << R*U << endl;

   double lambda1, lambda2, lambda3;
   Tensor2Gen EVec;
   U.eigen( &lambda1, &lambda2, &lambda3, &EVec );
   cout << "eigen values of U: "
        << lambda1 << " " << lambda2 << " " << lambda3 << endl << endl;

   Tensor2Sym F1( lambda1*selfDyadic(EVec(1)) +
                  lambda2*selfDyadic(EVec(2)) +
                  lambda3*selfDyadic(EVec(3)) );

   cout << "spectral composition (should be U)" << endl
        << setw(12) << F1 << endl;
   cout << "R*U*Rt (=V)" << endl << setw(12) << multASAt( R, U ) << endl;

   cout << "Ft*F (general mult)" << endl << setw(12) << F.transpose()*F << endl;
   cout << "C=Ft*F (multAtA)" << endl << setw(12) << multAtA(F) << endl;
   cout << "C=square(U)" << endl << setw(12) << square(U) << endl;

   Tensor2Sym V;
   polarLeft( F, &V, &R );

   cout << "polar decomposition F=VR" << endl;
   cout << "R" << endl << setw(12) << R << endl;
   cout << "V" << endl << setw(12) << V << endl;
   cout << "V.R" << endl << setw(12) << V*R << endl;
   cout << "Rt*V*R (=U)" << endl << setw(12) << multAtSA( R, V ) << endl;

   cout << "log(V)" << endl << setw(12) << log(V) << endl;

   // test article Franca
   Tensor2Sym Cs, U_inv; 

   Cs(1,1) = 400;
   Cs(2,2) = 0.9999000;
   Cs(3,3) = 0.0001999867;
   Cs(2,1) = 0;
   Cs(3,2) = 0.009998333;
   Cs(3,1) = 0;

   findroots( Cs, &U, &U_inv );  
   cout << setprecision(8);
   cout << "C" << endl << setw(14) << Cs << endl << "roots" << endl
        << "U" << endl << setw(14) << U << endl
        << "U_inv" << endl << setw(14) << U_inv << endl;
   cout << setprecision(5);

   // Rigid rotation -90 degrees:
   F(1,1) = 0; F(1,2) = -1; F(1,3) = 0;
   F(2,1) = 1; F(2,2) =  0; F(2,3) = 0;
   F(3,1) = 0; F(3,2) =  0; F(3,3) = 1;

   polarRight( F, &R, &U );

   cout << "F(Rigid rotation -90 degrees" << endl << setw(12) << F << endl;
   cout << "R" << endl << setw(12) << R << endl;
   cout << "U" << endl << setw(12) << U << endl;
   cout << "R.U" << endl << setw(12) << R*U << endl;

   cout << "log(U)" << endl << setw(12) << log(U) << endl;

   Tensor2Sym sig(0);
   sig(2,1) = 1.0;

   Tensor2Gen Q(0);
   Q(1,1) =  1/sqrt(2.0);
   Q(1,2) =  1/sqrt(2.0);
   Q(2,1) = -1/sqrt(2.0);
   Q(2,2) =  1/sqrt(2.0);

   cout << "sig(0)" << endl << setw(12) << sig << endl;
   cout << "rotated stress sig(45) = Q*sig(0)*Qt" << endl;
   cout << "sig(45)" << endl << setw(12) << multASAt( Q, sig ) << endl;

//  elasticity with general 2nd and 4th order tensors

   Tensor4Gen E4, II;
   Tensor2Gen sigma, eps, e;
   double  E, nu, G, K, p;

   E=200000.0; nu=0.3;

   G = E/(2*(1+nu));
   K = E/(3*(1-2*nu));

   Tensor4Gen H(1);

   Tensor2Gen I(1);
   II = dyadic(I,I);

   E4 = (2*G)*H + (K-2*G/3)*II;

   cout << "E4 = 2G H + (K-2G/3)II" << endl << setw(12) << E4 << endl;

   eps = 0;
   eps(1,1) = 0.001;
   eps(2,2) = -0.0003;
   eps(3,3) = -0.0003;

   eps(1,2) = 0.002;
   eps(2,1) = 0.002;
   eps(1,3) = 0.004;
   eps(2,3) = -0.015;
   eps(3,1) = 0.004;
   eps(3,2) = -0.015;

   cout << "eps" << endl << setw(12) << eps << endl;

   sigma = E4*eps;
   cout << "sigma 1 = E4*eps" << endl << setw(12) << sigma << endl;

   e = dev(eps);

   p = -K*trace(eps);
   sigma = (2*G)*e - p*I;
   cout << "sigma 2 = 2G*dev(eps) + K trace(eps)" << endl
        << setw(12) << sigma << endl;

// elasticity with symmetric 2nd and 4th order tensors

   Tensor4DSym E4s, IIs;
   Tensor2Sym sigmas, epss, es;

   Tensor4DSym Hs(1);
   Tensor2Sym  Is(1);

   IIs = selfDyadic(Is);

   E4s = (2*G)*Hs + (K-2*G/3)*IIs;

   cout << "E4 symmetric" << endl << setw(12) << E4s << endl;

   epss = 0;
   epss(1,1) = 0.001;
   epss(2,2) = -0.0003;
   epss(3,3) = -0.0003;

   epss(2,1) = 0.002;
   epss(3,1) = 0.004;
   epss(3,2) = -0.015;

   cout << "eps symmetric" << endl << setw(12) << epss << endl;

   sigmas = E4s*epss;
   cout << "sigma symmetric 1" << endl << setw(12) << sigmas << endl;

   es = dev(epss);

   p = -K*trace(epss);
   sigmas = (2*G)*es - p*Is;
   cout << "sigma symmetric 2" << endl << setw(12) << sigmas << endl;

   Tensor4LSym Hl(1);
   Tensor4LSym E4l;
   E4l = (2*G)*Hl + (K-2*G/3)*dyadic(Is,Is);
   sigmas = E4l*epss;
   cout << "sigma symmetric 3" << endl << setw(12) << sigmas << endl;

   cout << "double contractions" << endl;
   cout << "sigma:sigma   " << doubleContraction( sigma, sigma ) << endl;
   cout << "sigmas:sigmas " << doubleContraction( sigmas, sigmas ) << endl;
   cout << "sigma:sigmas  "
        << doubleContraction( (const Tensor2&)sigma, (const Tensor2&)sigmas )
        << endl << endl;

   E4s(2,1,1,1) = 140;
   E4s(2,1,2,2) = 150;
   E4s(2,1,3,3) = 160;
   E4s(3,2,1,1) = 170;
   E4s(3,2,2,2) = 145;
   E4s(3,2,3,3) = 155;
   E4s(3,2,2,1) = 165;
   E4s(3,1,1,1) = 175;
   E4s(3,1,2,2) = 185;
   E4s(3,1,3,3) = 130;
   E4s(3,1,2,1) = 140;
   E4s(3,1,3,2) = 150;

   cout << "filled symmetric 4th order tensor E4s" << endl
        << setw(12) << E4s << endl;
   cout << "symmetric mult" << endl << setw(12) << E4s*epss << endl;
   E4 = E4s;
   eps = epss;
   cout << "nonsymmetric mult E4*eps" << endl << setw(12) << E4*eps << endl;
   cout << "nonsymmetric mult eps*E4" << endl << setw(12) << eps*E4 << endl;

   a(1) = a(2) = 0.707107; a(3) = 0;
   b(1) = -0.707107; b(2) = 0.707107; b(3) = 0;
   cout << "vector product" << endl << setw(12) << a << endl
        << setw(12) << b << endl << setw(12) << vectorProduct(a,b) << endl;

   double I1, I2, I3;
   epss.invariants( &I1, &I2, &I3 );

   epss.eigen( &lambda1, &lambda2, &lambda3 );
   I1 = lambda1+lambda2+lambda3;
   I2 = lambda1*lambda2 + lambda2*lambda3 + lambda1*lambda3;
   I3 = lambda1*lambda2*lambda3;

   cout << "eigenvalues epss " << setw(12) << lambda1 << setw(12) << lambda2
                               << setw(12) << lambda3 << endl << endl;
   cout << "invariants epss            "
        << setw(12) << I1 << setw(12) << I2 << setw(12) << I3 << endl;
   cout << "invariants via eigenvalues "
        << setw(12) << I1 << setw(12) << I2 << setw(12) << I3 << endl << endl;

   I2 = 0.5*( pow(trace(epss),2) - trace( epss*epss ) );
   cout << "second invariant as: 0.5( trace(A)^2 - trace(A^2) ) "
        << setw(12) << I2 << endl << endl;


   cout << "check conversion to and from fortan style vector and matrix"
        << endl << endl;

   cout << "epss as tensor" << endl << setw(12) << epss << endl;
   double epsvec[6];
   epss.getVector( epsvec );
   cout << "epss as vector (shear * sqrt(2))" << endl;
   for ( int i=0; i<6; i++ ) cout << setw(12) << epsvec[i] << " ";
   cout << endl << endl;
   Tensor2Sym p2;
   p2.putVector( epsvec );
   cout << "and put as tensor again" << endl << setw(12) << p2 << endl;

   cout << "E4s as tensor" << endl << setw(12) << E4s << endl;
   double E4mat[36];
   E4s.getFortranMatrix( E4mat );
   cout << "E4s as matrix (shear/normal * sqrt(2), shear/shear * 2)" << endl;
   for ( int i=0; i<6; i++ )
   {
      for ( int j=0; j<6; j++ )
         cout << setw(12) << E4mat[i+6*j] << " ";
      cout << endl;
   }
   cout << endl;
   Tensor4DSym ps4;
   ps4.putFortranMatrix( E4mat );
   cout << "and put as tensor again" << endl << setw(12) << ps4 << endl;

   Tensor4LSym p4;
   p4.putFortranMatrix( E4mat );
   cout << "and put as Tensor4LSym" << endl << setw(12) << p4 << endl;
   double p4mat[36];
   p4.getFortranMatrix( p4mat );
   cout << "Tensor4LSym put as matrix (shear/normal * sqrt(2), shear/shear * 2)"
        << endl;
   for ( int i=0; i<6; i++ )
   {
      for ( int j=0; j<6; j++ )
         cout << setw(12) << p4mat[i+6*j] << " ";
      cout << endl;
   }
   cout << endl;

//   getchar();  // for use with ms-windows

   return 0;
}
