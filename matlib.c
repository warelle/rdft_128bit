#include <stdio.h>
#include "matlib.h"

void vec_add_float128(__float128 a[MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 c[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
      c[i] = a[i] + b[i];
}
void vec_sub_float128(__float128 a[MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 c[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
    c[i] = a[i] - b[i];
}
void vec_add_complex128(__complex128 a[MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 c[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
      c[i] = a[i] + b[i];
}
void vec_sub_complex128(__complex128 a[MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 c[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
    c[i] = a[i] - b[i];
}

void mat_vec_dot_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 c[MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    c[i] = 0.0Q;
    for(j=0; j<MATRIX_SIZE; j++){
      c[i] += a[i][j]*b[j];
    }
  }
}
void mat_vec_dot_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 c[MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    COMPLEX_ASSIGN(c[i], 0.0Q, 0.0Q);
    for(j=0; j<MATRIX_SIZE; j++){
      c[i] += a[i][j]*b[j];
    }
  }
}

void mat_add_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE], __float128 c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      c[i][j] = a[i][j] + b[i][j];
}
void mat_sub_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE], __float128 c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      c[i][j] = a[i][j] - b[i][j];
}
void mat_mul_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE][MATRIX_SIZE], __float128 c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j,k;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      c[i][j] = 0;
      for(k=0; k<MATRIX_SIZE; k++){
        c[i][j] += a[i][k]*b[k][j];
      }
    }
  }
}


void mat_add_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE], __complex128 c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      c[i][j] = a[i][j] + b[i][j];
}
void mat_sub_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE], __complex128 c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      c[i][j] = a[i][j] - b[i][j];
}
void mat_mul_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE][MATRIX_SIZE], __complex128 c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j,k;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      c[i][j] = 0;
      for(k=0; k<MATRIX_SIZE; k++){
        c[i][j] += a[i][k]*b[k][j];
      }
    }
  }
}

void print_vector_float128(__float128 f[MATRIX_SIZE]){
  int i;
  printf("[");
  for(i=0; i<MATRIX_SIZE; i++){
    print_float128(f[i]);
    printf(" ");
  }
  printf("]");
  printf("\n");
}
void print_vector_complex128(__complex128 f[MATRIX_SIZE]){
  int i;
  printf("[");
  for(i=0; i<MATRIX_SIZE; i++){
    print_complex128(f[i]);
    printf(" ");
  }
  printf("]");
  printf("\n");
}

void print_matrix_float128(__float128 f[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  printf("[");
  for(i=0; i<MATRIX_SIZE; i++){
    printf("[");
    for(j=0; j<MATRIX_SIZE; j++){
      print_float128(f[i][j]);
      printf(" ");
    }
    printf("]");
  }
  printf("]");
  printf("\n");
}
void print_matrix_complex128(__complex128 f[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  printf("[");
  for(i=0; i<MATRIX_SIZE; i++){
    printf("[");
    for(j=0; j<MATRIX_SIZE; j++){
      print_complex128(f[i][j]);
      printf(" ");
    }
    printf("]");
  }
  printf("]");
  printf("\n");
}
