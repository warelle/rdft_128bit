#include "gen.h"
#include "lib.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

// return: random value
double uniform(){
  static int init_flg = 0;
  if(!init_flg){
    init_flg = 1;
    srand((unsigned)time(NULL));
  }
  return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}
double rand_normal_double(double m, double s){
    double z = sqrt(-2.0*log(uniform()))*sin(2.0*M_PI*uniform());
    return m + s*z;
}
__float128 rand_normal_float_128(__float128 m, __float128 s){
    __float128 z = sqrtq(-2.0*logq(uniform()))*sinq(2.0*M_PIq*uniform());
    return m + s*z;
}
__float128 random_float_128(__float128 range){
  __float128 r = 0.0Q;
  r += uniform();
  r *= 1.0e-10;
  r += uniform();
  r *= 1.0e-10;
  r += uniform();
  r *= range;
  return uniform() > 0.5 ? r : -r;
}
__complex128 random_complex_128(__float128 range){
  __complex128 r;
  COMPLEX_ASSIGN(r, random_float_128(range),random_float_128(range));
  return r;
}

// return: vec
void generate_vector_float128(__float128 vec[MATRIX_SIZE], __float128 range){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
    vec[i] = random_float_128(range);
}
void generate_vector_complex128(__complex128 vec[MATRIX_SIZE], __float128 range){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
    vec[i] = random_complex_128(range);
}

// return: mat
void generate_matrix_float128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 range, int band_size){
//  gen_random_float128(mat,range);
//  gen_diag_big_float128(mat,range);
//  gen_arrowhead_float128(mat,range);
//  gen_identity_float128(mat);
//  gen_band_no_side_float128(mat, range, band_size);
  gen_band_float128(mat, range, band_size);
}
void generate_matrix_complex128(__complex128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 range){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      mat[i][j] = random_complex_128(range);
}

// return: mat,vec
void generate_linear_system_float_128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 vec[MATRIX_SIZE], __float128 range, int band_size){
  generate_matrix_float128(mat, range, band_size);
  generate_vector_float128(vec, range);
}
void generate_linear_system_complex_128(__complex128 mat[MATRIX_SIZE][MATRIX_SIZE], __complex128 vec[MATRIX_SIZE], __float128 range){
  generate_matrix_complex128(mat, range);
  generate_vector_complex128(vec, range);
}


/*--- specific type of matrices ---*/
void gen_identity_float128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      if(i == j)
        mat[i][j] = 1.0Q;
      else
        mat[i][j] = 0.0Q;
}
void gen_random_float128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 range){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      mat[i][j] = random_float_128(range);
}
void gen_arrowhead_float128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 range){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      if(i == j || i==0 || j==0)
        mat[i][j] = random_float_128(range);
      else
        mat[i][j] = 0.0Q;
}
void gen_diag_big_float128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 range){
  int i,j;
  __float128 big_f = 1000000000000000000.0Q;
  __float128 small_f = 1.0Q;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      if(i == j)
        mat[i][j] = random_float_128(range*big_f);
      else
        mat[i][j] = random_float_128(range*small_f);
}
// [NOTE] BAND_SIZE will be 2*band_size-1
void gen_band_no_side_float128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 range, int band_size){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    if(i<band_size){
      mat[i][0] = random_float_128(range);
      mat[0][i] = random_float_128(range);
    }else{
      mat[i][0] = 0.0Q;
      mat[0][i] = 0.0Q;
    }
  }
  for(i=1; i<MATRIX_SIZE; i++)
    for(j=1; j<MATRIX_SIZE; j++)
      if(mat[i-1][j-1] != 0.0Q)
        mat[i][j] = random_float_128(range);
      else
        mat[i][j] = 0.0Q;
}
void gen_band_float128(__float128 mat[MATRIX_SIZE][MATRIX_SIZE], __float128 range, int band_size){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    if(i<band_size || (MATRIX_SIZE-band_size) < i){
      mat[i][0] = random_float_128(range);
      mat[0][i] = random_float_128(range);
    }else{
      mat[i][0] = 0.0Q;
      mat[0][i] = 0.0Q;
    }
  }
  for(i=1; i<MATRIX_SIZE; i++)
    for(j=1; j<MATRIX_SIZE; j++)
      if(mat[i-1][j-1] != 0.0Q)
        mat[i][j] = random_float_128(range);
      else
        mat[i][j] = 0.0Q;
}
