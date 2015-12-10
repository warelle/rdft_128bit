#ifndef MYRDFT_H
#define MYRDFT_H

#include <quadmath.h>
#include <complex>
#include "test.h"

// 64 bit
void dft_matrix_complex_double(std::complex<double> f[MATRIX_SIZE][MATRIX_SIZE]);
void gauss_matrix_double(double f[MATRIX_SIZE][MATRIX_SIZE]);
void r_matrix_complex_double(std::complex<double> r[MATRIX_SIZE][MATRIX_SIZE]);

// 128 bit
void dft_matrix_complex_128(__complex128 f[MATRIX_SIZE][MATRIX_SIZE]);
void dht_matrix(__float128 f[MATRIX_SIZE][MATRIX_SIZE]);
void gauss_matrix_float128(__float128 f[MATRIX_SIZE][MATRIX_SIZE]);
void r_matrix_complex_128(__complex128 r[MATRIX_SIZE][MATRIX_SIZE]);
void r_real_matrix(__float128 r[MATRIX_SIZE][MATRIX_SIZE]);

#endif
