#ifndef MYGIVENS_H
#define MYGIVENS_H

#include "test.h"
#include "matlib.h"
#include <complex>

//
//   | i  j
// ---------
// i | c -s
// j | s  c
//
typedef struct{
  int i;
  int j;
  std::complex<double> c;
  std::complex<double> s;
} givens_matrix;

typedef struct givens_matrix_list{
  givens_matrix *gm;
  struct givens_matrix_list *next;
}givens_matrix_list;

givens_matrix* create_givens_matrix(int i, int j, double theta);
void delete_givens_matrix(givens_matrix *gm);

givens_matrix_list *r_givens_matrix_double();

givens_matrix_list *create_givens_matrix_list();
void delete_givens_matrix_list(givens_matrix_list *gml);

void convert_givens_matrix_complex_double(givens_matrix *gm, std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE]);

void ins_r_givens_matrix_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE]);

void mat_vec_dot_givens_complex_double(givens_matrix_list *gml, std::complex<double> a[MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE]);
//void mat_mul_givens_left_complex_double(givens_matrix_list *gml, std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE]);
void mat_mul_givens_right_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], givens_matrix_list *gml, std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE]);

#endif
