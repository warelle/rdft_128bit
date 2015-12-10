#include "eigen.h"

#include <Eigen/Core>
#include <Eigen/SVD>

using namespace Eigen;
using namespace std;

complex<double> ac[MATRIX_SIZE][MATRIX_SIZE];

__complex128 down_cast_complex_128(__complex128 a){
  __complex128 r;
  __float128 f1,f2;
  double a1 = REALPART(a), a2 = IMAGPART(a);
  complex<double> rr(a1, a2);
  f1 = real(rr);
  f2 = imag(rr);
  COMPLEX_ASSIGN(r, f1,f2);
  return r;
}
void down_cast_vec_complex_128(__complex128 a[MATRIX_SIZE], __complex128 b[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++){
      __complex128 tmp = a[i];
      b[i] = down_cast_complex_128(tmp);
  }
}
void down_cast_mat_complex_128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      __complex128 tmp = a[i][j];
      b[i][j] = down_cast_complex_128(tmp);
    }
  }
}

void down_cast_complex128_to_complex(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], complex<double> b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      double a1 = REALPART(a[i][j]), a2 = IMAGPART(a[i][j]);
      complex<double> rr(a1, a2);
      b[i][j] = rr;
    }
  }
}

// Jacobi SVD for complex maybe has some bug
//double condition_number_complex_128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], int k, double *sfrak){
//  int i,j;
//  down_cast_complex128_to_complex(a, ac);
//
//  MatrixXcd mat;
//  mat.resize(k,k);
//
//  for(i=0; i<k; i++)
//    for(j=0; j<k; j++)
//      mat(i,j) = ac[i][j];
//
//  JacobiSVD<MatrixXcd> svd(mat);
//
//  auto sval = svd.singularValues();
//  if(sfrak != NULL)
//    *sfrak = sval[sval.size()-1];
//
//
//  return sval[0]/sval[sval.size()-1];
//}

double condition_number(double a[MATRIX_SIZE][MATRIX_SIZE], double *sa){
  MatrixXd mat =  Map<Matrix<double, Dynamic, Dynamic> >(&(a[0][0]), MATRIX_SIZE, MATRIX_SIZE);
  JacobiSVD<MatrixXd> svd(mat);

  auto sval = svd.singularValues();
  if(sa != NULL)
    *sa = sval[sval.size()-1];

  return sval[0]/sval[sval.size()-1];
}
