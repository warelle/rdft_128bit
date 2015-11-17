#include <stdio.h>
#include "lib.h"

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

void cast_mat_float128_to_double(__float128 a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      b[i][j] = (double)a[i][j];
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
