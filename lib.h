#ifndef MYLIB_H
#define MYLIB_H

#include <quadmath.h>
#include <complex>
#include "test.h"

#define REALPART(z) (__real__(z))
#define IMAGPART(z) (__imag__(z))
#define COMPLEX_ASSIGN(z_, r_, i_) {__real__(z_) = (r_); __imag__(z_) = (i_);}

// norm
double vector_norm_double(double a[MATRIX_SIZE]);
__float128 vector_norm_float128(__float128 a[MATRIX_SIZE]);
__float128 vector_norm_complex128(__complex128 a[MATRIX_SIZE]);

// down accuracy
__float128 down_cast_float128(__float128 a);
void down_cast_vec_float128(__float128 a[MATRIX_SIZE], __float128 b[MATRIX_SIZE]);
void down_cast_mat_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE]);

// cast 64-128 or 128-64 float only
void cast_vec_float128_to_double(__float128 a[MATRIX_SIZE], double b[MATRIX_SIZE]);
void cast_mat_float128_to_double(__float128 a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE]);

// cast 64-128 or 128-64 complex only

// cast 64-128 or 128-64 mixed
void cast_vec_double_to_complex_double(double a[MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE]);
void cast_vec_complex_double_to_double(std::complex<double> a[MATRIX_SIZE], double b[MATRIX_SIZE]);
void cast_mat_double_to_complex_double(double a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE]);
void cast_mat_complex_double_to_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE]);

// cast 128 bit only
void cast_vec_float128_to_complex128(__float128 a[MATRIX_SIZE], __complex128 b[MATRIX_SIZE]);
void cast_vec_complex128_to_float128(__complex128 a[MATRIX_SIZE], __float128 b[MATRIX_SIZE]);
void cast_mat_float128_to_complex128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE]);
void cast_mat_complex128_to_float128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE]);

// print
void print_double(double d);
void print_float128(__float128 f);
void print_complex128(__complex128 d);

#endif
