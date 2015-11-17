#ifndef MYGEN_H
#define MYGEN_H

#include "test.h"

// return: random value
double uniform();
__float128 rand_normal(__float128 m, __float128 s);
__float128 random_float_128(__float128 range);
__complex128 random_complex_128(__float128 range);

// return: vec
void generate_vector_float128(__float128 vec[MATRIX_SIZE], __float128 range);
void generate_vector_complex128(__complex128 vec[MATRIX_SIZE], __float128 range);

// return: mat
void generate_matrix_float128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 range);
void generate_matrix_complex128(__complex128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 range);

// return: mat,vec
void generate_linear_system_float_128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 vec[MATRIX_SIZE], __float128 range);
void generate_linear_system_complex_128(__complex128 mat[MATRIX_SIZE][MATRIX_SIZE], __complex128 vec[MATRIX_SIZE], __float128 range);


#endif
