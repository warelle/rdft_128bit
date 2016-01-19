#include <cstdio>
#include <cmath>
#include <complex>
#include "test.h"
#include "lib.h"
#include "lu.h"
#include "gen.h"
#include "rdft.h"
#include "matlib.h"
#include "solve.h"

#include "eigen.h"

#define NONE                       0
#define RDFT                       1
#define RDFT_ITERATION             2
#define DHT                        4
#define DHT_ITERATION              8
#define GAUSS                      16
#define GAUSS_ITERATION            32
#define PP                         64
#define PP_ITERATION               128
#define RDFT_PERM                  256
#define RDFT_PERM_ITERATION        512
#define RDFT_GIVENS                1024
#define RDFT_GIVENS_ITERATION      2048
#define RDFT_GIVENS_TWO            4096
#define RDFT_GIVENS_TWO_ITERATION  8192
#define RDFT_BOTH_GIVENS           16384
#define RDFT_BOTH_GIVENS_ITERATION 32768
#define ALL                       131071


__float128 a[MATRIX_SIZE][MATRIX_SIZE];
__float128 b[MATRIX_SIZE];
__float128 x[MATRIX_SIZE];

__complex128 aa[MATRIX_SIZE][MATRIX_SIZE];
__complex128 bb[MATRIX_SIZE];
__complex128 xx[MATRIX_SIZE];
__complex128 xi[MATRIX_SIZE];
__complex128 xa[MATRIX_SIZE];

// for calc condition number
double double_a[MATRIX_SIZE][MATRIX_SIZE];

// double
std::complex<double> dc_aa[MATRIX_SIZE][MATRIX_SIZE];
std::complex<double> dc_bb[MATRIX_SIZE];
std::complex<double> dc_xx[MATRIX_SIZE];
std::complex<double> dc_xi[MATRIX_SIZE];
std::complex<double> dc_xa[MATRIX_SIZE];

double d_a[MATRIX_SIZE][MATRIX_SIZE];
double d_b[MATRIX_SIZE];
double d_x[MATRIX_SIZE];
double d_x_rdft[MATRIX_SIZE];
double d_x_rdft_iter[MATRIX_SIZE];
double d_x_rdft_iter_another[MATRIX_SIZE];
double d_x_rdft_perm[MATRIX_SIZE];
double d_x_rdft_perm_iter[MATRIX_SIZE];
double d_x_rdft_perm_iter_another[MATRIX_SIZE];
double d_x_rdft_givens[MATRIX_SIZE];
double d_x_rdft_givens_iter[MATRIX_SIZE];
double d_x_rdft_givens_iter_another[MATRIX_SIZE];
double d_x_rdft_givens_two[MATRIX_SIZE];
double d_x_rdft_givens_two_iter[MATRIX_SIZE];
double d_x_rdft_givens_two_iter_another[MATRIX_SIZE];
double d_x_rdft_both_givens[MATRIX_SIZE];
double d_x_rdft_both_givens_iter[MATRIX_SIZE];
double d_x_rdft_both_givens_iter_another[MATRIX_SIZE];
double d_x_rdht[MATRIX_SIZE];
double d_x_rdht_iter[MATRIX_SIZE];
double d_x_rdht_iter_another[MATRIX_SIZE];
double d_x_gauss[MATRIX_SIZE];
double d_x_gauss_iter[MATRIX_SIZE];
double d_x_gauss_iter_another[MATRIX_SIZE];

double d_x_rdft_perm_dif[MATRIX_SIZE];
double d_x_rdft_givens_dif[MATRIX_SIZE];
double d_x_rdft_givens_two_dif[MATRIX_SIZE];
double d_x_rdft_both_givens_dif[MATRIX_SIZE];
double d_x_rdht_dif[MATRIX_SIZE];
double d_x_gauss_dif[MATRIX_SIZE];

int d_rdft_err;
int d_rdft_perm_err;
int d_rdft_givens_err;
int d_rdft_givens_two_err;
int d_rdft_both_givens_err;
int d_gauss_err;

void run_64bit(int dat, int opt, int exe, int band_size, int x_axis){
  generate_linear_system_float_128(a,x,100.0Q, band_size);
  cast_mat_float128_to_double(a,d_a);
  cast_vec_float128_to_double(x,d_x);

  down_cast_mat_float128(a,a);
  down_cast_vec_float128(x,x);
  mat_vec_dot_float128(a,x,b);

  cast_vec_float128_to_double(b,d_b);

  cast_mat_double_to_complex_double(d_a,dc_aa);
  cast_vec_double_to_complex_double(d_b,dc_bb);

  if(exe & (RDFT_PERM | RDFT_PERM_ITERATION)){
    solve_with_rdft_perm_iteration_complex_double(dc_aa, dc_bb, dc_xx, dc_xi, dc_xa, &d_rdft_perm_err);
  }
  if(exe & (RDFT_GIVENS | RDFT_GIVENS_ITERATION)){
    solve_with_rdft_givens_iteration_complex_double(dc_aa, dc_bb, dc_xx, dc_xi, dc_xa, &d_rdft_givens_err);
  }
  if(exe & (RDFT_GIVENS_TWO | RDFT_GIVENS_TWO_ITERATION)){
    solve_with_rdft_givens_two_iteration_complex_double(dc_aa, dc_bb, dc_xx, dc_xi, dc_xa, &d_rdft_givens_two_err);
  }
  if(exe & (RDFT_BOTH_GIVENS | RDFT_BOTH_GIVENS_ITERATION)){
    solve_with_rdft_both_givens_iteration_complex_double(dc_aa, dc_bb, dc_xx, dc_xi, dc_xa, &d_rdft_both_givens_err);
  }
  if(exe & (GAUSS | GAUSS_ITERATION)){
    solve_with_gauss_iteration_double(d_a, d_b, d_x_gauss, d_x_gauss_iter, d_x_gauss_iter_another, &d_gauss_err);
  }

  cast_mat_float128_to_double(a,double_a);

  if(opt == 0){
    printf("condition number:");
    printf("%f\n", condition_number(double_a, NULL));
    printf("RDFTM          :");
    print_int(d_rdft_perm_err);
    printf("\n");
    printf("RDFTG          :");
    print_int(d_rdft_givens_err);
    printf("\n");
    printf("RDFTG          :");
    print_int(d_rdft_givens_two_err);
    printf("\n");
    printf("RDFTBG         :");
    print_int(d_rdft_both_givens_err);
    printf("\n");
    printf("GAUSS          :");
    print_int(d_gauss_err);
    printf("\n");
  }else if(opt == 1){
    printf("%d ", x_axis);
    printf("%f ", condition_number(double_a,NULL));
    print_int(d_rdft_perm_err);
    printf(" ");
    print_int(d_rdft_givens_err);
    printf(" ");
    print_int(d_rdft_givens_two_err);
    printf(" ");
    print_int(d_rdft_both_givens_err);
    printf(" ");
    print_int(d_gauss_err);
    printf("\n");
  }else if(opt == 2){
    printf("condition number:");
    printf("%f\n", condition_number(double_a,NULL));
    if(exe & RDFT_PERM){
      printf("RDFT PERM          :");
      print_int(d_rdft_perm_err);
      printf("\n");
    }
    if(exe & RDFT_GIVENS){
      printf("RDFT GVN           :");
      print_int(d_rdft_givens_err);
      printf("\n");
    }
    if(exe & RDFT_GIVENS_TWO){
      printf("RDFT GTWO          :");
      print_int(d_rdft_givens_two_err);
      printf("\n");
    }
    if(exe & RDFT_BOTH_GIVENS){
      printf("RDFT BOTH          :");
      print_int(d_rdft_both_givens_err);
      printf("\n");
    }
    if(exe & GAUSS){
      printf("GAUSS              :");
      print_int(d_gauss_err);
      printf("\n");
    }
  }
}

int main(){
  int i,j;
  //int exe = RDFT_PERM | RDFT_PERM_ITERATION | RDFT_GIVENS | RDFT_GIVENS_ITERATION | RDFT_GIVENS_TWO | RDFT_GIVENS_TWO_ITERATION | RDFT_BOTH_GIVENS | RDFT_BOTH_GIVENS_ITERATION | GAUSS | GAUSS_ITERATION | PP | PP_ITERATION | LIB;
  //int exe = RDFT | RDFT_ITERATION | RDFT_PERM | RDFT_PERM_ITERATION | RDFT_GIVENS | RDFT_GIVENS_ITERATION | GAUSS | GAUSS_ITERATION | PP | PP_ITERATION;
  int exe = ALL;
  //for(i=0; i<MATRIX_SIZE; i++){
  for(i=0; i<100; i++){
    //for(j=0; j<4; j++)
      run_64bit(i,1,exe, i+1, MATRIX_SIZE);
    fprintf(stderr, "%d ", i);
    //fprintf(stderr, "%d ", MATRIX_SIZE);
  }

  return 0;
}
