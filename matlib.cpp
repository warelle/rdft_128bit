#include <stdio.h>
#include "matlib.h"

// ----- vector atith ----- //
void vec_add_double(double a[MATRIX_SIZE], double b[MATRIX_SIZE], double c[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
      c[i] = a[i] + b[i];
}
void vec_sub_double(double a[MATRIX_SIZE], double b[MATRIX_SIZE], double c[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
    c[i] = a[i] - b[i];
}
void vec_add_complex_double(std::complex<double> a[MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
      c[i] = a[i] + b[i];
}
void vec_sub_complex_double(std::complex<double> a[MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
    c[i] = a[i] - b[i];
}
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

// ----- matrix vector atith ----- //
void mat_vec_dot_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double c[MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    c[i] = 0.0;
    for(j=0; j<MATRIX_SIZE; j++){
      c[i] += a[i][j]*b[j];
    }
  }
}
void mat_vec_dot_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    c[i] = std::complex<double>(0,0);
    for(j=0; j<MATRIX_SIZE; j++){
      c[i] += a[i][j]*b[j];
    }
  }
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

// ----- matrix atith ----- //
void mat_add_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE], double c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      c[i][j] = a[i][j] + b[i][j];
}
void mat_sub_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE], double c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      c[i][j] = a[i][j] - b[i][j];
}
void mat_mul_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE], double c[MATRIX_SIZE][MATRIX_SIZE]){
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
void mat_add_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      c[i][j] = a[i][j] + b[i][j];
}
void mat_sub_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      c[i][j] = a[i][j] + b[i][j];
}
void mat_mul_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> c[MATRIX_SIZE][MATRIX_SIZE]){
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

// ----- matrix transpose ----- //
void mat_transpose_double(double a[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
        double tmp = a[i][j];
        a[i][j] = a[j][i];
        a[j][i] = tmp;
    }
  }
}
void mat_transpose_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
        std::complex<double> tmp = a[i][j];
        a[i][j] = a[j][i];
        a[j][i] = tmp;
    }
  }
}
void mat_transpose_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
        __float128 tmp = a[i][j];
        a[i][j] = a[j][i];
        a[j][i] = tmp;
    }
  }
}
void mat_transpose_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE]){
   int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
        __complex128 tmp = a[i][j];
        a[i][j] = a[j][i];
        a[j][i] = tmp;
    }
  }
}

// ----- vector matrix print ----- //
void print_vector_double(double f[MATRIX_SIZE]){
  int i;
  printf("[");
  for(i=0; i<MATRIX_SIZE; i++){
    printf("%f",f[i]);
    printf(" ");
  }
  printf("]");
  printf("\n");
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

void print_matrix_double(double f[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  printf("[");
  for(i=0; i<MATRIX_SIZE; i++){
    printf("[");
    for(j=0; j<MATRIX_SIZE; j++){
      printf("%f",f[i][j]);
      printf(" ");
    }
    printf("]");
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
