#ifndef MYMATLIB_H
#define MYMATLIB_H

#include <quadmath.h>
#include "lib.h"
#include "test.h"

// print vector, matrix
void print_vector_double(double f[MATRIX_SIZE]);
void print_vector_complex_double(std::complex<double> f[MATRIX_SIZE]);
void print_vector_float128(__float128 f[MATRIX_SIZE]);
void print_vector_complex128(__complex128 f[MATRIX_SIZE]);
void print_matrix_double(double f[MATRIX_SIZE][MATRIX_SIZE]);
void print_matrix_complex_double(std::complex<double> f[MATRIX_SIZE][MATRIX_SIZE]);
void print_matrix_float128(__float128 f[MATRIX_SIZE][MATRIX_SIZE]);
void print_matrix_complex128(__complex128 f[MATRIX_SIZE][MATRIX_SIZE]);

// c= a+-b
// input a,b
// return c
void vec_add_double(double a[MATRIX_SIZE], double b[MATRIX_SIZE], double c[MATRIX_SIZE]);
void vec_sub_double(double a[MATRIX_SIZE], double b[MATRIX_SIZE], double c[MATRIX_SIZE]);
void vec_add_complex_double(std::complex<double> a[MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE]);
void vec_sub_complex_double(std::complex<double> a[MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE]);
void vec_add_float128(__float128 a[MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 c[MATRIX_SIZE]);
void vec_sub_float128(__float128 a[MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 c[MATRIX_SIZE]);
void vec_add_complex128(__complex128 a[MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 c[MATRIX_SIZE]);
void vec_sub_complex128(__complex128 a[MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 c[MATRIX_SIZE]);


// c= A*b
// input a,b
// return c
void mat_vec_dot_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double c[MATRIX_SIZE]);
void mat_vec_dot_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE]);
void mat_vec_dot_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 c[MATRIX_SIZE]);
void mat_vec_dot_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 c[MATRIX_SIZE]);


// C = A+B, A-B, A*B
// input a,b
// return c
void mat_add_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE], double c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_sub_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE], double c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_mul_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE], double c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_add_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_sub_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_mul_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_add_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE], __float128 c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_sub_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE], __float128 c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_mul_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE], __float128 c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_add_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE], __complex128 c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_sub_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE], __complex128 c[MATRIX_SIZE][MATRIX_SIZE]);
void mat_mul_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE], __complex128 c[MATRIX_SIZE][MATRIX_SIZE]);

// C = C^T
void mat_transpose_double(double a[MATRIX_SIZE][MATRIX_SIZE]);
void mat_transpose_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE]);
void mat_transpose_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE]);
void mat_transpose_complex128(__float128 a[MATRIX_SIZE][MATRIX_SIZE]);


#endif
