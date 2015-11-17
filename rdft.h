#ifndef MYRDFT_H
#define MYRDFT_H

#include <quadmath.h>
#include "test.h"

void dft_matrix(__complex128 f[MATRIX_SIZE][MATRIX_SIZE]);
void dht_matrix(__float128 f[MATRIX_SIZE][MATRIX_SIZE]);
void gauss_matrix(__float128 f[MATRIX_SIZE][MATRIX_SIZE]);
void r_matrix(__complex128 r[MATRIX_SIZE][MATRIX_SIZE]);
void r_real_matrix(__float128 r[MATRIX_SIZE][MATRIX_SIZE]);

#endif
