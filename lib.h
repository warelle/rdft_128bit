#ifndef MYLIB_H
#define MYLIB_H

#include <quadmath.h>
#include "test.h"

#define REALPART(z) (__real__(z))
#define IMAGPART(z) (__imag__(z))
#define COMPLEX_ASSIGN(z_, r_, i_) {__real__(z_) = (r_); __imag__(z_) = (i_);}

__float128 vector_norm_float128(__float128 a[MATRIX_SIZE]);
__float128 vector_norm_complex128(__complex128 a[MATRIX_SIZE]);

void cast_vec_float128_to_complex128(__float128 a[MATRIX_SIZE], __complex128 b[MATRIX_SIZE]);
void cast_vec_complex128_to_float128(__complex128 a[MATRIX_SIZE], __float128 b[MATRIX_SIZE]);
void cast_mat_float128_to_double(__float128 a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE]);
void cast_mat_float128_to_complex128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE]);
void cast_mat_complex128_to_float128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE]);

void print_float128(__float128 f);
void print_complex128(__complex128 d);

#endif
