#ifndef MYLU_H
#define MYLU_H

#include <math.h>
#include <complex>
#include "lib.h"
#include "test.h"

// LU decomposition
//  input : m
//  return: l
//  return: u
void lu_double(double m[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE]);
void lu_complex_double(std::complex<double> m[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE]);
void lu_float128(__float128 m[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE]);
void lu_complex128(__complex128 m[MATRIX_SIZE][MATRIX_SIZE], __complex128 l[MATRIX_SIZE][MATRIX_SIZE], __complex128 u[MATRIX_SIZE][MATRIX_SIZE]);

// LU decomposition with partial pivoting
//  input : m
//  return: l
//  return: u
void lu_partial_pivot_double(double m[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], int p[MATRIX_SIZE]);
void lu_partial_pivot_float128(__float128 m[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE], int p[MATRIX_SIZE]);

// LU solver
//  input : l and b,  u and y
//  return: y      ,  x
void l_step_double(double l[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double y[MATRIX_SIZE]);
void u_step_double(double u[MATRIX_SIZE][MATRIX_SIZE], double y[MATRIX_SIZE], double x[MATRIX_SIZE]);
void l_step_complex_double(std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> y[MATRIX_SIZE]);
void u_step_complex_double(std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> y[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]);
void l_step_float128(__float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 y[MATRIX_SIZE]);
void u_step_float128(__float128 u[MATRIX_SIZE][MATRIX_SIZE], __float128 y[MATRIX_SIZE], __float128 x[MATRIX_SIZE]);
void l_step_complex128(__complex128 l[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 y[MATRIX_SIZE]);
void u_step_complex128(__complex128 u[MATRIX_SIZE][MATRIX_SIZE], __complex128 y[MATRIX_SIZE], __complex128 x[MATRIX_SIZE]);

#endif
