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
__float128 x_rdft_iter_another[MATRIX_SIZE];
__float128 x_rdht[MATRIX_SIZE];
__float128 x_rdht_iter[MATRIX_SIZE];
__float128 x_rdht_iter_another[MATRIX_SIZE];
__float128 x_gauss[MATRIX_SIZE];
__float128 x_gauss_iter[MATRIX_SIZE];
__float128 x_gauss_iter_another[MATRIX_SIZE];
__float128 x_pp[MATRIX_SIZE];
__float128 x_pp_iter[MATRIX_SIZE];
__float128 x_pp_iter_another[MATRIX_SIZE];

__complex128 aa[MATRIX_SIZE][MATRIX_SIZE];
__complex128 bb[MATRIX_SIZE];
__complex128 xx[MATRIX_SIZE];
__complex128 xi[MATRIX_SIZE];
__complex128 xa[MATRIX_SIZE];

__float128 x_rdft_dif[MATRIX_SIZE];
__float128 x_rdft_iter_dif[MATRIX_SIZE];
__float128 x_rdft_iter_another_dif[MATRIX_SIZE];
__float128 x_rdht_dif[MATRIX_SIZE];
__float128 x_rdht_iter_dif[MATRIX_SIZE];
__float128 x_rdht_iter_another_dif[MATRIX_SIZE];
__float128 x_gauss_dif[MATRIX_SIZE];
__float128 x_gauss_iter_dif[MATRIX_SIZE];
__float128 x_gauss_iter_another_dif[MATRIX_SIZE];
__float128 x_pp_dif[MATRIX_SIZE];
__float128 x_pp_iter_dif[MATRIX_SIZE];
__float128 x_pp_iter_another_dif[MATRIX_SIZE];

__float128 rdft_err;
__float128 rdft_iter_err;
__float128 rdft_iter_another_err;
__float128 rdht_err;
__float128 rdht_iter_err;
__float128 rdht_iter_another_err;
__float128 gauss_err;
__float128 gauss_iter_err;
__float128 gauss_iter_another_err;
__float128 pp_err;
__float128 pp_iter_err;
__float128 pp_iter_another_err;

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
double d_x_rdht[MATRIX_SIZE];
double d_x_rdht_iter[MATRIX_SIZE];
double d_x_rdht_iter_another[MATRIX_SIZE];
double d_x_gauss[MATRIX_SIZE];
double d_x_gauss_iter[MATRIX_SIZE];
double d_x_gauss_iter_another[MATRIX_SIZE];
double d_x_pp[MATRIX_SIZE];
double d_x_pp_iter[MATRIX_SIZE];
double d_x_pp_iter_another[MATRIX_SIZE];

double d_x_rdft_dif[MATRIX_SIZE];
double d_x_rdft_iter_dif[MATRIX_SIZE];
double d_x_rdft_iter_another_dif[MATRIX_SIZE];
double d_x_rdht_dif[MATRIX_SIZE];
double d_x_rdht_iter_dif[MATRIX_SIZE];
double d_x_rdht_iter_another_dif[MATRIX_SIZE];
double d_x_gauss_dif[MATRIX_SIZE];
double d_x_gauss_iter_dif[MATRIX_SIZE];
double d_x_gauss_iter_another_dif[MATRIX_SIZE];
double d_x_pp_dif[MATRIX_SIZE];
double d_x_pp_iter_dif[MATRIX_SIZE];
double d_x_pp_iter_another_dif[MATRIX_SIZE];

double d_rdft_err;
double d_rdft_iter_err;
double d_rdft_iter_another_err;
double d_rdht_err;
double d_rdht_iter_err;
double d_rdht_iter_another_err;
double d_gauss_err;
double d_gauss_iter_err;
double d_gauss_iter_another_err;
double d_pp_err;
double d_pp_iter_err;
double d_pp_iter_another_err;

void run_64bit(int dat, int opt, int exe){
  generate_linear_system_float_128(a,x,100.0Q);
  cast_mat_float128_to_double(a,d_a);
  cast_vec_float128_to_double(x,d_x);

  down_cast_mat_float128(a,a);
  down_cast_vec_float128(x,x);
  mat_vec_dot_float128(a,x,b);

  cast_vec_float128_to_double(b,d_b);

  cast_mat_double_to_complex_double(d_a,dc_aa);
  cast_vec_double_to_complex_double(d_b,dc_bb);

  if(exe & (RDFT | RDFT_ITERATION)){
    solve_with_rdft_iteration_complex_double(dc_aa, dc_bb, dc_xx, dc_xi, dc_xa);
    cast_vec_complex_double_to_double(dc_xx, d_x_rdft);
    cast_vec_complex_double_to_double(dc_xi, d_x_rdft_iter);
    cast_vec_complex_double_to_double(dc_xa, d_x_rdft_iter_another);
    vec_sub_double(d_x,d_x_rdft,d_x_rdft_dif);
    vec_sub_double(d_x,d_x_rdft_iter,d_x_rdft_iter_dif);
    vec_sub_double(d_x,d_x_rdft_iter_another,d_x_rdft_iter_another_dif);
    d_rdft_err = vector_norm_double(d_x_rdft_dif);
    d_rdft_iter_err = vector_norm_double(d_x_rdft_iter_dif);
    d_rdft_iter_another_err = vector_norm_double(d_x_rdft_iter_another_dif);
  }
  if(exe & (GAUSS | GAUSS_ITERATION)){
    solve_with_gauss_iteration_double(d_a, d_b, d_x_gauss, d_x_gauss_iter, d_x_gauss_iter_another);
    vec_sub_double(d_x,d_x_gauss,d_x_gauss_dif);
    vec_sub_double(d_x,d_x_gauss_iter,d_x_gauss_iter_dif);
    vec_sub_double(d_x,d_x_gauss_iter_another,d_x_gauss_iter_another_dif);
    d_gauss_err = vector_norm_double(d_x_gauss_dif);
    d_gauss_iter_err = vector_norm_double(d_x_gauss_iter_dif);
    d_gauss_iter_another_err = vector_norm_double(d_x_gauss_iter_another_dif);
  }
  if(exe & (PP | PP_ITERATION)){
    solve_with_partial_pivot_double(d_a,d_b,d_x_pp, d_x_pp_iter, d_x_pp_iter_another);
    vec_sub_double(d_x,d_x_pp,d_x_pp_dif);
    vec_sub_double(d_x,d_x_pp_iter,d_x_pp_iter_dif);
    vec_sub_double(d_x,d_x_pp_iter_another,d_x_pp_iter_another_dif);
    d_pp_err   = vector_norm_double(d_x_pp_dif);
    d_pp_iter_err   = vector_norm_double(d_x_pp_iter_dif);
    d_pp_iter_another_err   = vector_norm_double(d_x_pp_iter_another_dif);
  }

  cast_mat_float128_to_double(a,double_a);

  if(opt == 0){
    printf("condition number:");
    printf("%f\n", condition_number(double_a, NULL));
    printf("RDFT           :");
    print_double(d_rdft_err);
    printf("\n");
    printf("RDFT iteration :");
    print_double(d_rdft_iter_err);
    printf("\n");
    printf("RDFT iteration :");
    print_double(d_rdft_iter_another_err);
    printf("\n");
    printf("GAUSS          :");
    print_double(d_gauss_err);
    printf("\n");
    printf("GAUSS iteration:");
    print_double(d_gauss_iter_err);
    printf("\n");
    printf("GAUSS iteration:");
    print_double(d_gauss_iter_another_err);
    printf("\n");
    printf("Partial Pivot  :");
    print_double(d_pp_err);
    printf("\n");
    printf("Partial Pivot  :");
    print_double(d_pp_iter_err);
    printf("\n");
    printf("Partial Pivot  :");
    print_double(d_pp_iter_another_err);
    printf("\n");
  }else if(opt == 1){ // graph data
    if(condition_number(double_a,NULL) > 15000)
      return;
    printf("%f ", condition_number(double_a,NULL));
    print_double(d_rdft_err);
    printf(" ");
    print_double(d_rdft_iter_err);
    printf(" ");
    print_double(d_rdft_iter_another_err);
    printf(" ");
    print_double(d_gauss_err);
    printf(" ");
    print_double(d_gauss_iter_err);
    printf(" ");
    print_double(d_gauss_iter_another_err);
    printf(" ");
    print_double(d_pp_err);
    printf(" ");
    print_double(d_pp_iter_err);
    printf(" ");
    print_double(d_pp_iter_another_err);
    printf("\n");
  }else if(opt == 2){
    printf("condition number:");
    printf("%f\n", condition_number(double_a,NULL));
    if(exe & RDFT){
      printf("RDFT           :");
      print_double(d_rdft_err);
      printf("\n");
    }
    if(exe & RDFT_ITERATION){
      printf("RDFT iteration :");
      print_double(d_rdft_iter_err);
      printf("\n");
      printf("RDFT iteration :");
      print_double(d_rdft_iter_another_err);
      printf("\n");
    }
    if(exe & GAUSS){
      printf("GAUSS          :");
      print_double(d_gauss_err);
      printf("\n");
    }
    if(exe & GAUSS_ITERATION){
      printf("GAUSS iteration:");
      print_double(d_gauss_iter_err);
      printf("\n");
      printf("GAUSS iteration:");
      print_double(d_gauss_iter_another_err);
      printf("\n");
    }
    if(exe & PP){
      printf("Partial Pivot  :");
      print_double(d_pp_err);
      printf("\n");
    }
    if(exe & PP_ITERATION){
      printf("Partial Pivot  :");
      print_double(d_pp_iter_err);
      printf("\n");
      printf("Partial Pivot  :");
      print_double(d_pp_iter_another_err);
      printf("\n");
    }
  }
}

void run_128bit(int dat, int opt, int exe){
  generate_linear_system_float_128(a,x,100.0Q);
  mat_vec_dot_float128(a,x,b);

  // [CARE] this must change the result
  //down_cast_vec_float_128(b,b);

  cast_mat_float128_to_complex128(a,aa);
  cast_vec_float128_to_complex128(b,bb);

  if(exe & (RDFT | RDFT_ITERATION)){
    solve_with_rdft_iteration_complex128(aa, bb, xx, xi, xa);
    cast_vec_complex128_to_float128(xx, x_rdft);
    cast_vec_complex128_to_float128(xi, x_rdft_iter);
    cast_vec_complex128_to_float128(xa, x_rdft_iter_another);
    vec_sub_float128(x,x_rdft,x_rdft_dif);
    vec_sub_float128(x,x_rdft_iter,x_rdft_iter_dif);
    vec_sub_float128(x,x_rdft_iter_another,x_rdft_iter_another_dif);
    rdft_err = vector_norm_float128(x_rdft_dif);
    rdft_iter_err = vector_norm_float128(x_rdft_iter_dif);
    rdft_iter_another_err = vector_norm_float128(x_rdft_iter_another_dif);
  }
  if(exe & (DHT | DHT_ITERATION)){
    solve_with_rdht_iteration_float128(a, b, x_rdht, x_rdht_iter, x_rdht_iter_another);
    vec_sub_float128(x,x_rdht,x_rdht_dif);
    vec_sub_float128(x,x_rdht_iter,x_rdht_iter_dif);
    vec_sub_float128(x,x_rdht_iter_another,x_rdht_iter_another_dif);
    rdht_err = vector_norm_float128(x_rdht_dif);
    rdht_iter_err = vector_norm_float128(x_rdht_iter_dif);
    rdht_iter_another_err = vector_norm_float128(x_rdht_iter_another_dif);
  }
  if(exe & (GAUSS | GAUSS_ITERATION)){
    solve_with_gauss_iteration_float128(a, b, x_gauss, x_gauss_iter, x_gauss_iter_another);
    vec_sub_float128(x,x_gauss,x_gauss_dif);
    vec_sub_float128(x,x_gauss_iter,x_gauss_iter_dif);
    vec_sub_float128(x,x_gauss_iter_another,x_gauss_iter_another_dif);
    gauss_err = vector_norm_float128(x_gauss_dif);
    gauss_iter_err = vector_norm_float128(x_gauss_iter_dif);
    gauss_iter_another_err = vector_norm_float128(x_gauss_iter_another_dif);
  }
  if(exe & (PP | PP_ITERATION)){
    solve_with_partial_pivot_float128(a,b,x_pp, x_pp_iter, x_pp_iter_another);
    vec_sub_float128(x,x_pp,x_pp_dif);
    vec_sub_float128(x,x_pp_iter,x_pp_iter_dif);
    vec_sub_float128(x,x_pp_iter_another,x_pp_iter_another_dif);
    pp_err   = vector_norm_float128(x_pp_dif);
    pp_iter_err   = vector_norm_float128(x_pp_iter_dif);
    pp_iter_another_err   = vector_norm_float128(x_pp_iter_another_dif);
  }

  cast_mat_float128_to_double(a,double_a);

  if(opt == 0){
    printf("condition number:");
    printf("%f\n", condition_number(double_a, NULL));
    printf("RDFT           :");
    print_float128(rdft_err);
    printf("\n");
    printf("RDFT iteration :");
    print_float128(rdft_iter_err);
    printf("\n");
    printf("RDFT iteration :");
    print_float128(rdft_iter_another_err);
    printf("\n");
    printf("RDHT           :");
    print_float128(rdht_err);
    printf("\n");
    printf("RDHT iteration :");
    print_float128(rdht_iter_err);
    printf("\n");
    printf("RDHT iteration :");
    print_float128(rdht_iter_another_err);
    printf("\n");
    printf("GAUSS          :");
    print_float128(gauss_err);
    printf("\n");
    printf("GAUSS iteration:");
    print_float128(gauss_iter_err);
    printf("\n");
    printf("GAUSS iteration:");
    print_float128(gauss_iter_another_err);
    printf("\n");
    printf("Partial Pivot  :");
    print_float128(pp_err);
    printf("\n");
    printf("Partial Pivot  :");
    print_float128(pp_iter_err);
    printf("\n");
    printf("Partial Pivot  :");
    print_float128(pp_iter_another_err);
    printf("\n");
  }else if(opt == 1){ // graph data
    printf("%f ", condition_number(double_a,NULL));
    print_float128(rdft_err);
    printf(" ");
    print_float128(rdft_iter_err);
    printf(" ");
    print_float128(rdft_iter_another_err);
    printf(" ");
    print_float128(rdht_err);
    printf(" ");
    print_float128(rdht_iter_err);
    printf(" ");
    print_float128(rdht_iter_another_err);
    printf(" ");
    print_float128(gauss_err);
    printf(" ");
    print_float128(gauss_iter_err);
    printf(" ");
    print_float128(gauss_iter_another_err);
    printf(" ");
    print_float128(pp_err);
    printf(" ");
    print_float128(pp_iter_err);
    printf(" ");
    print_float128(pp_iter_another_err);
    printf("\n");
  }else if(opt == 2){
    printf("condition number:");
    printf("%f\n", condition_number(double_a,NULL));
    if(exe & RDFT){
      printf("RDFT           :");
      print_float128(rdft_err);
      printf("\n");
    }
    if(exe & RDFT_ITERATION){
      printf("RDFT iteration :");
      print_float128(rdft_iter_err);
      printf("\n");
      printf("RDFT iteration :");
      print_float128(rdft_iter_another_err);
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
      printf("DHT iteration  :");
      print_float128(rdht_iter_another_err);
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
      printf("GAUSS iteration:");
      print_float128(gauss_iter_another_err);
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
      printf("Partial Pivot  :");
      print_float128(pp_iter_another_err);
      printf("\n");
    }
  }
}

extern __complex128 c_f[MATRIX_SIZE][MATRIX_SIZE];
extern __complex128 c_r[MATRIX_SIZE][MATRIX_SIZE];
extern __complex128 c_fr[MATRIX_SIZE][MATRIX_SIZE];
extern __complex128 c_fra[MATRIX_SIZE][MATRIX_SIZE];

//void singular_value_test(){
//  int i;
//  double cond_a, cond_fra, sfrak=0, sa=0, minsfrak=0, maxcondfrak=0;
//
//  generate_linear_system_float_128(a,x,100.0Q);
//
//  cast_mat_float128_to_complex128(a,aa);
//  cast_mat_float128_to_double(a,double_a);
//
//  dft_matrix(c_f);
//  r_matrix(c_r);
//  mat_mul_complex128(c_f,c_r,c_fr);
//  mat_mul_complex128(c_fr,aa,c_fra);
//
//  cond_a   = condition_number(double_a, &sa);
//  printf("cond(FRA_k),max,  s_n(FRA_k), min\n");
//  for(i=1; i<=MATRIX_SIZE; i++){
//    cond_fra = condition_number_complex_128(c_fra, i, &sfrak);
//    if(cond_fra > maxcondfrak || i == 1)
//      maxcondfrak = cond_fra;
//    if(sfrak < minsfrak || i == 1)
//      minsfrak = sfrak;
//    printf("%f,%f ,%f,%f \n", cond_fra, maxcondfrak, sfrak, minsfrak);
//  }
//
//  printf("s_n(A)     :%f\n", sa);
//  printf("s_n(FRA_k) :%f\n", minsfrak);
//  printf("cond(A)    :%f\n", cond_a);
//  printf("cond(FRA_k):%f\n", maxcondfrak);
//  printf("%f <= %f \n", maxcondfrak/cond_a, sa/minsfrak);
//}

int main(){
  int i;
  for(i=0; i<100; i++){
    int exe = RDFT | RDFT_ITERATION | DHT | DHT_ITERATION | GAUSS | GAUSS_ITERATION | PP | PP_ITERATION;
    run_64bit(i,1,exe);
  }

  //for(i=0; i<50; i++){
  //  singular_value_test();
  //}

  return 0;
}
