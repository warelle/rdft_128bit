#include "rdft.h"
#include "lib.h"
#include "matlib.h"
#include "gen.h"

#include <random>
#include <vector>

int for_perm[MATRIX_SIZE];
int init_for_perm_flg = 0;

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

void init_for_perm(){
  int i;
  for(i=0; i<MATRIX_SIZE; i++){
    for_perm[i] = i;
  }
  init_for_perm_flg = 1;
}

void r_perm_matrix_complex_double(std::complex<double> r[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  double nn = MATRIX_SIZE;
  std::random_device rd;
	std::mt19937 mt(rd());

  if(! init_for_perm_flg){
    init_for_perm();
  }

  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      r[i][j] = std::complex<double>(0,0);
    }
  }

  std::vector<int> perm(for_perm,for_perm+MATRIX_SIZE);

  for(i=0; i<MATRIX_SIZE; i++){
    double rval = uniform()*M_PI;
    int idx = mt() % (MATRIX_SIZE - i);
    if(rval < 0.0)
      rval = -rval;
    r[i][perm[idx]] = std::complex<double>( cos(-2.0*M_PI*rval/nn), sin(-2.0*M_PI*rval/nn) );
    perm.erase(perm.begin()+idx);
  }
}
void r_perm_matrix_complex_128(__complex128 r[MATRIX_SIZE][MATRIX_SIZE]){
   int i,j;
  __float128 nn = MATRIX_SIZE;
  std::random_device rd;
	std::mt19937 mt(rd());

  if(! init_for_perm_flg){
    init_for_perm();
  }

  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      COMPLEX_ASSIGN(r[i][j], 0.0Q, 0.0Q);
    }
  }

  std::vector<int> perm(for_perm,for_perm+MATRIX_SIZE);

  for(i=0; i<MATRIX_SIZE; i++){
    __float128 rval = random_float_128(M_PIq);
    int idx = mt() % (MATRIX_SIZE - i);
     if(rval < 0.0Q)
       rval = -rval;
     COMPLEX_ASSIGN(r[i][perm[idx]], cosq(-2.0Q*M_PIq*rval/nn), sinq(-2.0Q*M_PIq*rval/nn));
     perm.erase(perm.begin()+idx);
  }
}

void r_real_matrix(__float128 r[MATRIX_SIZE][MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++){
    r[i][i] = 1.0Q;
  }
}
