#include <cstdio>
#include "lib.h"


// ----- norm ----- //
double vector_norm_double(double a[MATRIX_SIZE]){
   int i=0;
   double r = 0.0;
   for(i=0; i<MATRIX_SIZE; i++){
      r += a[i]*a[i];
    }
    r = sqrt(r);
    return r;
}
__float128 vector_norm_float128(__float128 a[MATRIX_SIZE]){
   int i=0;
   __float128 r = 0.0Q;
   for(i=0; i<MATRIX_SIZE; i++){
      r += a[i]*a[i];
    }
    r = sqrtq(r);
    return r;
}
__float128 vector_norm_complex128(__complex128 a[MATRIX_SIZE]){
  int i;
  __float128 r = 0.0Q;
  for(i=0; i<MATRIX_SIZE; i++){
    r += REALPART(a[i]*conjq(a[i]));
  }
  r = sqrtq(r);
  return r;
}

// ----- down accuracy ----- //
__float128 down_cast_float128(__float128 a){
  __float128 r;
  double rr = (double)a;
  r = (__float128)rr;
  return r;
}
void down_cast_vec_float128(__float128 a[MATRIX_SIZE], __float128 b[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++){
      double tmp = a[i];
      b[i] = (__float128)tmp;
  }
}
void down_cast_mat_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      double tmp = a[i][j];
      b[i][j] = (__float128)tmp;
    }
  }
}


// ----- cast 64-128 or 128-64 float only ----- //
void cast_vec_float128_to_double(__float128 a[MATRIX_SIZE], double b[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++){
    b[i] = (double)a[i];
  }
}
void cast_mat_float128_to_double(__float128 a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      b[i][j] = (double)a[i][j];
    }
  }
}

// ----- cast 64-128 or 128-64 complex only ----- //

// ----- cast 64-128 or 128-64 mixed ----- //
void cast_vec_double_to_complex_double(double a[MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
      b[i] = std::complex<double>(a[i],0);
}
void cast_mat_double_to_complex_double(double a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      b[i][j] = std::complex<double>(a[i][j], 0);
}
void cast_vec_complex_double_to_double(std::complex<double> a[MATRIX_SIZE], double b[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
      b[i] = std::real(a[i]);
}
void cast_mat_complex_double_to_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      b[i][j] = std::real(a[i][j]);
}


// ----- cast 128bit only ----- //
void cast_vec_float128_to_complex128(__float128 a[MATRIX_SIZE], __complex128 b[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
      COMPLEX_ASSIGN(b[i], a[i], 0.0Q);
}
void cast_vec_complex128_to_float128(__complex128 a[MATRIX_SIZE], __float128 b[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
      b[i] = REALPART(a[i]);
}

void cast_mat_float128_to_complex128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      COMPLEX_ASSIGN(b[i][j], a[i][j], 0.0Q);
}
void cast_mat_complex128_to_float128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      b[i][j] = REALPART(a[i][j]);
}


// ----- count zero ----- //
int count_zero_mat_double(double a[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j,r;
  r=0;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      if(-1.0e-30 <= a[i][j] && a[i][j] <= 1.0e-30)
        r++;
  return r;
}
int count_zero_mat_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j,r;
  r=0;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      if(-1.0e-30 <= std::norm(a[i][j]) && std::norm(a[i][j]) <= 1.0e-30)
        r++;
  return r;
}



// ----- print ----- //
void print_int(int d){
  printf("%d",d);
}
void print_double(double d){
  printf("%.25f",d);
}
void print_float128(__float128 f){
  char c[100];
  quadmath_snprintf(c,80,"%.40Qf", f);
  printf("%s",c);
}

void print_complex128(__complex128 d){
  char c[100];
  quadmath_snprintf(c,80,"%.40Qf", REALPART(d));
  printf("%s",c);
  quadmath_snprintf(c,80,"%.40Qf", IMAGPART(d));
  if(c[0] != '-')
    printf("+");
  printf("%sj",c);
}
