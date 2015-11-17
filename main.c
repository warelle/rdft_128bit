#include <stdio.h>
#include <math.h>
#include "test.h"
#include "lib.h"
#include "lu.h"
#include "gen.h"
#include "matlib.h"
#include "solve.h"

#include "eigen.h"

#define NONE               0
#define RDFT               1
#define RDFT_ITERATION     2
#define DHT                4
#define DHT_ITERATION      8
#define GAUSS              16
#define GAUSS_ITERATION    32
#define PP                 64
#define PP_ITERATION       128


__float128 a[MATRIX_SIZE][MATRIX_SIZE];
__float128 b[MATRIX_SIZE];
__float128 x[MATRIX_SIZE];
__float128 x_rdft[MATRIX_SIZE];
__float128 x_rdft_iter[MATRIX_SIZE];
__float128 x_rdht[MATRIX_SIZE];
__float128 x_rdht_iter[MATRIX_SIZE];
__float128 x_gauss[MATRIX_SIZE];
__float128 x_gauss_iter[MATRIX_SIZE];
__float128 x_pp[MATRIX_SIZE];
__float128 x_pp_iter[MATRIX_SIZE];

__complex128 aa[MATRIX_SIZE][MATRIX_SIZE];
__complex128 bb[MATRIX_SIZE];
__complex128 xx[MATRIX_SIZE];
__complex128 xi[MATRIX_SIZE];

__float128 x_rdft_dif[MATRIX_SIZE];
__float128 x_rdft_iter_dif[MATRIX_SIZE];
__float128 x_rdht_dif[MATRIX_SIZE];
__float128 x_rdht_iter_dif[MATRIX_SIZE];
__float128 x_gauss_dif[MATRIX_SIZE];
__float128 x_gauss_iter_dif[MATRIX_SIZE];
__float128 x_pp_dif[MATRIX_SIZE];
__float128 x_pp_iter_dif[MATRIX_SIZE];

__float128 rdft_err;
__float128 rdft_iter_err;
__float128 rdht_err;
__float128 rdht_iter_err;
__float128 gauss_err;
__float128 gauss_iter_err;
__float128 pp_err;
__float128 pp_iter_err;

double double_a[MATRIX_SIZE][MATRIX_SIZE];

void run(int dat, int opt, int exe){
  generate_linear_system_float_128(a,x,100.0Q);
  mat_vec_dot_float128(a,x,b);

  cast_mat_float128_to_complex128(a,aa);
  cast_vec_float128_to_complex128(b,bb);

  if(exe & (RDFT | RDFT_ITERATION)){
    solve_with_rdft_iteration(aa, bb, xx, xi);
    cast_vec_complex128_to_float128(xx, x_rdft);
    cast_vec_complex128_to_float128(xi, x_rdft_iter);
    vec_sub_float128(x,x_rdft,x_rdft_dif);
    vec_sub_float128(x,x_rdft_iter,x_rdft_iter_dif);
    rdft_err = vector_norm_float128(x_rdft_dif);
    rdft_iter_err = vector_norm_float128(x_rdft_iter_dif);
  }
  if(exe & (DHT | DHT_ITERATION)){
    solve_with_rdht_iteration(a, b, x_rdht, x_rdht_iter);
    vec_sub_float128(x,x_rdht,x_rdht_dif);
    vec_sub_float128(x,x_rdht_iter,x_rdht_iter_dif);
    rdht_err = vector_norm_float128(x_rdht_dif);
    rdht_iter_err = vector_norm_float128(x_rdht_iter_dif);
  }
  if(exe & (GAUSS | GAUSS_ITERATION)){
    solve_with_gauss_iteration(a, b, x_gauss, x_gauss_iter);
    vec_sub_float128(x,x_gauss,x_gauss_dif);
    vec_sub_float128(x,x_gauss_iter,x_gauss_iter_dif);
    gauss_err = vector_norm_float128(x_gauss_dif);
    gauss_iter_err = vector_norm_float128(x_gauss_iter_dif);
  }
  if(exe & (PP | PP_ITERATION)){
    solve_with_partial_pivot(a,b,x_pp, x_pp_iter);
    vec_sub_float128(x,x_pp,x_pp_dif);
    vec_sub_float128(x,x_pp_iter,x_pp_iter_dif);
    pp_err   = vector_norm_float128(x_pp_dif);
    pp_iter_err   = vector_norm_float128(x_pp_iter_dif);
  }

  cast_mat_float128_to_double(a,double_a);

  if(opt == 0){
    printf("condition number:");
    printf("%f\n", condition_number(double_a));
    printf("RDFT           :");
    print_float128(rdft_err);
    printf("\n");
    printf("RDFT iteration :");
    print_float128(rdft_iter_err);
    printf("\n");
    printf("RDHT           :");
    print_float128(rdht_err);
    printf("\n");
    printf("RDHT iteration :");
    print_float128(rdht_iter_err);
    printf("\n");
    printf("GAUSS          :");
    print_float128(gauss_err);
    printf("\n");
    printf("GAUSS iteration:");
    print_float128(gauss_iter_err);
    printf("\n");
    printf("Partial Pivot  :");
    print_float128(pp_err);
    printf("\n");
    printf("Partial Pivot  :");
    print_float128(pp_iter_err);
    printf("\n");
  }else if(opt == 1){ // graph data
    printf("%f ", condition_number(double_a));
    print_float128(rdft_err);
    printf(" ");
    print_float128(rdft_iter_err);
    printf(" ");
    print_float128(rdht_err);
    printf(" ");
    print_float128(rdht_iter_err);
    printf(" ");
    print_float128(gauss_err);
    printf(" ");
    print_float128(gauss_iter_err);
    printf(" ");
    print_float128(pp_err);
    printf(" ");
    print_float128(pp_iter_err);
    printf("\n");
  }else if(opt == 2){
    printf("condition number:");
    printf("%f\n", condition_number(double_a));
    if(exe & RDFT){
      printf("RDFT           :");
      print_float128(rdft_err);
      printf("\n");
    }
    if(exe & RDFT_ITERATION){
      printf("RDFT iteration :");
      print_float128(rdft_iter_err);
      printf("\n");
    }
    if(exe & DHT){
      printf("DHT            :");
      print_float128(rdht_err);
      printf("\n");
    }
    if(exe & DHT_ITERATION){
      printf("DHT iteration  :");
      print_float128(rdht_iter_err);
      printf("\n");
    }
    if(exe & GAUSS){
      printf("GAUSS          :");
      print_float128(gauss_err);
      printf("\n");
    }
    if(exe & GAUSS_ITERATION){
      printf("GAUSS iteration:");
      print_float128(gauss_iter_err);
      printf("\n");
    }
    if(exe & PP){
      printf("Partial Pivot  :");
      print_float128(pp_err);
      printf("\n");
    }
    if(exe & PP_ITERATION){
      printf("Partial Pivot  :");
      print_float128(pp_iter_err);
      printf("\n");
    }
  }
}

int main(){
  int i;
  for(i=0; i<100; i++){
    int exe = RDFT | RDFT_ITERATION | DHT | DHT_ITERATION | GAUSS | GAUSS_ITERATION | PP | PP_ITERATION;
    run(i,1,exe);
    printf("\n");
  }

  return 0;
}
