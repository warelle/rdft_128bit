#include "lu.h"

double d_mat[MATRIX_SIZE][MATRIX_SIZE];
std::complex<double> cd_mat[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_mat[MATRIX_SIZE][MATRIX_SIZE];
__complex128 c_mat[MATRIX_SIZE][MATRIX_SIZE];

void lu_double(double m[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j,k;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      d_mat[i][j] = m[i][j];

  for(i=0; i<MATRIX_SIZE-1; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
      d_mat[j][i] /= d_mat[i][i];
      for(k=i+1; k<MATRIX_SIZE; k++){
        d_mat[j][k] -= d_mat[j][i]*d_mat[i][k];
      }
    }
  }

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++){
    l[i][i] = 1.0Q;
    for(j=0; j<MATRIX_SIZE; j++)
      if(i>j)
        l[i][j] = d_mat[i][j];
      else
        u[i][j] = d_mat[i][j];
  }
}
void lu_complex_double(std::complex<double> m[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j,k;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      cd_mat[i][j] = m[i][j];

  for(i=0; i<MATRIX_SIZE-1; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
      cd_mat[j][i] /= cd_mat[i][i];
      for(k=i+1; k<MATRIX_SIZE; k++){
        cd_mat[j][k] -= cd_mat[j][i]*cd_mat[i][k];
      }
    }
  }

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++){
    l[i][i] = 1.0Q;
    for(j=0; j<MATRIX_SIZE; j++)
      if(i>j)
        l[i][j] = cd_mat[i][j];
      else
        u[i][j] = cd_mat[i][j];
  }
}
void lu_float128(__float128 m[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j,k;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      f_mat[i][j] = m[i][j];

  for(i=0; i<MATRIX_SIZE-1; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
      f_mat[j][i] /= f_mat[i][i];
      for(k=i+1; k<MATRIX_SIZE; k++){
        f_mat[j][k] -= f_mat[j][i]*f_mat[i][k];
      }
    }
  }

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++){
    l[i][i] = 1.0Q;
    for(j=0; j<MATRIX_SIZE; j++)
      if(i>j)
        l[i][j] = f_mat[i][j];
      else
        u[i][j] = f_mat[i][j];
  }
}

void lu_complex128(__complex128 m[MATRIX_SIZE][MATRIX_SIZE], __complex128 l[MATRIX_SIZE][MATRIX_SIZE], __complex128 u[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j,k;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      c_mat[i][j] = m[i][j];

  for(i=0; i<MATRIX_SIZE-1; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
      c_mat[j][i] /= c_mat[i][i];
      for(k=i+1; k<MATRIX_SIZE; k++){
        c_mat[j][k] -= c_mat[j][i]*c_mat[i][k];
      }
    }
  }

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++){
    COMPLEX_ASSIGN(l[i][i], 1.0Q, 0.0Q);
    for(j=0; j<MATRIX_SIZE; j++)
      if(i>j)
        l[i][j] = c_mat[i][j];
      else
        u[i][j] = c_mat[i][j];
  }
}

void swap_double(double m[MATRIX_SIZE][MATRIX_SIZE], int i, int j){
    double tmp[MATRIX_SIZE];
    int k;
    for(k=0; k<MATRIX_SIZE; k++){
      tmp[k] = m[i][k];
      m[i][k] = m[j][k];
      m[j][k] = tmp[k];
    }
}
int pivot_double(double m[MATRIX_SIZE][MATRIX_SIZE], int k){
  int i;
  int piv = -1;
  int line = k;
  for(i=k; i<MATRIX_SIZE; i++){
    if(fabs(m[i][k]) > piv){
      piv = fabs(m[i][k]);
      line = i;
    }
  }
  if(line != k){
    swap_double(m, line, k);
  }
  return line;
}
void lu_partial_pivot_double(double m[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], int p[MATRIX_SIZE]){
  int i,j,k;
  double mat[MATRIX_SIZE][MATRIX_SIZE];

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      mat[i][j] = m[i][j];

  for(i=0; i<MATRIX_SIZE-1; i++){
    p[i] = pivot_double(mat, i);
    for(j=i+1; j<MATRIX_SIZE; j++){
      mat[j][i] /= mat[i][i];
      for(k=i+1; k<MATRIX_SIZE; k++){
        mat[j][k] -= mat[j][i]*mat[i][k];
      }
    }
  }

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++){
    l[i][i] = 1.0Q;
    for(j=0; j<MATRIX_SIZE; j++)
      if(i>j)
        l[i][j] = mat[i][j];
      else
        u[i][j] = mat[i][j];
  }
}

void swap_float128(__float128 m[MATRIX_SIZE][MATRIX_SIZE], int i, int j){
    __float128 tmp[MATRIX_SIZE];
    int k;
    for(k=0; k<MATRIX_SIZE; k++){
      tmp[k] = m[i][k];
      m[i][k] = m[j][k];
      m[j][k] = tmp[k];
    }
}
int pivot_float128(__float128 m[MATRIX_SIZE][MATRIX_SIZE], int k){
  int i;
  int piv = -1;
  int line = k;
  for(i=k; i<MATRIX_SIZE; i++){
    if(fabsq(m[i][k]) > piv){
      piv = fabsq(m[i][k]);
      line = i;
    }
  }
  if(line != k){
    swap_float128(m, line, k);
  }
  return line;
}
void lu_partial_pivot_float128(__float128 m[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE], int p[MATRIX_SIZE]){
  int i,j,k;
  __float128 mat[MATRIX_SIZE][MATRIX_SIZE];

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      mat[i][j] = m[i][j];

  for(i=0; i<MATRIX_SIZE-1; i++){
    p[i] = pivot_float128(mat, i);
    for(j=i+1; j<MATRIX_SIZE; j++){
      mat[j][i] /= mat[i][i];
      for(k=i+1; k<MATRIX_SIZE; k++){
        mat[j][k] -= mat[j][i]*mat[i][k];
      }
    }
  }

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++){
    l[i][i] = 1.0Q;
    for(j=0; j<MATRIX_SIZE; j++)
      if(i>j)
        l[i][j] = mat[i][j];
      else
        u[i][j] = mat[i][j];
  }
}

void l_step_double(double l[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double y[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      y[i] = b[i];

  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<i; j++){
      y[i] -= l[i][j]*y[j];
    }
    y[i] /= l[i][i];
  }
}
void u_step_double(double u[MATRIX_SIZE][MATRIX_SIZE], double y[MATRIX_SIZE], double x[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      x[i] = y[i];

  for(i=MATRIX_SIZE-1; i>=0; i--){
    for(j=MATRIX_SIZE-1; j>=i+1; j--){
      x[i] -= u[i][j]*x[j];
    }
    x[i] /= u[i][i];
  }
}

void l_step_complex_double(std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> y[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      y[i] = b[i];

  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<i; j++){
      y[i] -= l[i][j]*y[j];
    }
    y[i] /= l[i][i];
  }
}
void u_step_complex_double(std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> y[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      x[i] = y[i];

  for(i=MATRIX_SIZE-1; i>=0; i--){
    for(j=MATRIX_SIZE-1; j>=i+1; j--){
      x[i] -= u[i][j]*x[j];
    }
    x[i] /= u[i][i];
  }
}

void l_step_float128(__float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 y[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      y[i] = b[i];

  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<i; j++){
      y[i] -= l[i][j]*y[j];
    }
    y[i] /= l[i][i];
  }
}
void u_step_float128(__float128 u[MATRIX_SIZE][MATRIX_SIZE], __float128 y[MATRIX_SIZE], __float128 x[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      x[i] = y[i];

  for(i=MATRIX_SIZE-1; i>=0; i--){
    for(j=MATRIX_SIZE-1; j>=i+1; j--){
      x[i] -= u[i][j]*x[j];
    }
    x[i] /= u[i][i];
  }
}


void l_step_complex128(__complex128 l[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 y[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      y[i] = b[i];

  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<i; j++){
      y[i] -= l[i][j]*y[j];
    }
    y[i] /= l[i][i];
  }
}
void u_step_complex128(__complex128 u[MATRIX_SIZE][MATRIX_SIZE], __complex128 y[MATRIX_SIZE], __complex128 x[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      x[i] = y[i];

  for(i=MATRIX_SIZE-1; i>=0; i--){
    for(j=MATRIX_SIZE-1; j>=i+1; j--){
      x[i] -= u[i][j]*x[j];
    }
    x[i] /= u[i][i];
  }
}

