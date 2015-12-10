#include "rdft.h"
#include "lib.h"
#include "gen.h"

void dft_matrix_complex_double(std::complex<double> f[MATRIX_SIZE][MATRIX_SIZE]){
  int i, j;
  double nn = MATRIX_SIZE;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      f[i][j] = std::complex<double>(cosq(-2.0*M_PI*i*j/nn), sin(-2.0*M_PI*i*j/nn));
    }
  }
}
void dft_matrix_complex_128(__complex128 f[MATRIX_SIZE][MATRIX_SIZE]){
  int i, j;
  __float128 nn = MATRIX_SIZE;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      __float128 ii = i;
      __float128 jj = j;
      COMPLEX_ASSIGN(f[i][j], cosq(-2.0Q*M_PIq*ii*jj/nn), sinq(-2.0Q*M_PIq*ii*jj/nn));
    }
  }
}

void dht_matrix(__float128 f[MATRIX_SIZE][MATRIX_SIZE]){
  int i, j;
  __float128 nn = MATRIX_SIZE;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      __float128 ii = i;
      __float128 jj = j;
      f[i][j] = cosq(2.0Q*M_PIq*ii*jj/nn) + sinq(2.0Q*M_PIq*ii*jj/nn);
    }
  }
}

void gauss_matrix_double(double f[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      f[i][j] = rand_normal_double(0.0, 1.0);
    }
  }
}
void gauss_matrix_float128(__float128 f[MATRIX_SIZE][MATRIX_SIZE]){
  int i, j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      f[i][j] = rand_normal_float_128(0.0Q, 1.0Q);
    }
  }
}

void r_matrix_complex_double(std::complex<double> r[MATRIX_SIZE][MATRIX_SIZE]){
  int i;
  double nn = MATRIX_SIZE;
  for(i=0; i<MATRIX_SIZE; i++){
    double rval = uniform()*M_PI;
    if(rval < 0.0)
      rval = -rval;
    r[i][i] = std::complex<double>( cos(-2.0*M_PI*rval/nn), sin(-2.0*M_PI*rval/nn) );
  }
}
void r_matrix_complex_128(__complex128 r[MATRIX_SIZE][MATRIX_SIZE]){
   int i;
  __float128 nn = MATRIX_SIZE;
  for(i=0; i<MATRIX_SIZE; i++){
    __float128 rval = random_float_128(M_PIq);
     if(rval < 0.0Q)
       rval = -rval;
     COMPLEX_ASSIGN(r[i][i], cosq(-2.0Q*M_PIq*rval/nn), sinq(-2.0Q*M_PIq*rval/nn));
  }
}

void r_real_matrix(__float128 r[MATRIX_SIZE][MATRIX_SIZE]){
   int i;
  for(i=0; i<MATRIX_SIZE; i++){
    r[i][i] = 1.0Q;
  }
}
