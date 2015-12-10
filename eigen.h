#ifndef MYEIGEN_H
#define MYEIGEN_H

#include "test.h"
#include "lib.h"
#include <complex>

__complex128 down_cast_complex_128(__complex128 a);
void down_cast_vec_complex_128(__complex128 a[MATRIX_SIZE], __complex128 b[MATRIX_SIZE]);
void down_cast_mat_complex_128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE]);
void down_cast_complex128_to_complex(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE]);

double condition_number_complex_128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], int k, double *sfrak);
double condition_number(double a[MATRIX_SIZE][MATRIX_SIZE], double *sa);

#endif
