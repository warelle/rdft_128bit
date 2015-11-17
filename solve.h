#ifndef MYSOLVE_H
#define MYSOLVE_H

#include <quadmath.h>
#include "test.h"

// solve Ax = b
// input : a,b
// return: x
// [note] im afraid of neccesity of many information (r,fr,fra,...etc)
void solve_with_rdft_iteration(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE], __complex128 xi[MATRIX_SIZE]);

void solve_with_rdht_iteration(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE]);
void solve_with_gauss_iteration(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE]);

void solve_with_partial_pivot(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE]);

void iteration_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 l[MATRIX_SIZE][MATRIX_SIZE], __complex128 u[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE]);
void iteration_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE]);

#endif
