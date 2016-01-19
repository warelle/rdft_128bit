#ifndef MYSOLVE_H
#define MYSOLVE_H

#include <quadmath.h>
#include <complex>
#include "test.h"
#include "givens.h"

// solve Ax = b
// input : a,b
// return: x
// [note] im afraid of neccesity of many information (r,fr,fra,...etc)

// 64 bit
void solve_with_rdft_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE]);
void solve_with_rdft_perm_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter=NULL);
void solve_with_rdft_givens_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter=NULL);
void solve_with_rdft_givens_two_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter=NULL);
void solve_with_rdft_both_givens_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter=NULL);
void solve_with_gauss_iteration_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE], double xi[MATRIX_SIZE], double xia[MATRIX_SIZE], int *counter=NULL);
void solve_with_partial_pivot_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE], double xi[MATRIX_SIZE], double xia[MATRIX_SIZE]);

// 128 bit
void solve_with_rdft_iteration_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE], __complex128 xi[MATRIX_SIZE], __complex128 xia[MATRIX_SIZE]);
void solve_with_rdft_perm_iteration_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE], __complex128 xi[MATRIX_SIZE], __complex128 xia[MATRIX_SIZE]);
void solve_with_rdht_iteration_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE], __float128 xia[MATRIX_SIZE]);
void solve_with_gauss_iteration_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE], __float128 xia[MATRIX_SIZE]);
void solve_with_partial_pivot_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE], __float128 xia[MATRIX_SIZE]);

void iteration_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 l[MATRIX_SIZE][MATRIX_SIZE], __complex128 u[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE]);
void iteration_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE]);
void iteration_double(double a[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE]);
void iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]);

void iteration_complex128_another(__complex128 f[MATRIX_SIZE][MATRIX_SIZE], __complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 l[MATRIX_SIZE][MATRIX_SIZE], __complex128 u[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE]);
void iteration_float128_another(__float128 f[MATRIX_SIZE][MATRIX_SIZE], __float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE]);
void iteration_double_another(double f[MATRIX_SIZE][MATRIX_SIZE], double a[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE]);
void iteration_complex_double_another(std::complex<double> f[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]);

#endif
