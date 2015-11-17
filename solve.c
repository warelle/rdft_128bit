#include "solve.h"
#include "lu.h"
#include "rdft.h"
#include "matlib.h"
#include <stdio.h>

// solving
__complex128 c_f[MATRIX_SIZE][MATRIX_SIZE];
__complex128 c_r[MATRIX_SIZE][MATRIX_SIZE];
__complex128 c_fr[MATRIX_SIZE][MATRIX_SIZE];
__complex128 c_fra[MATRIX_SIZE][MATRIX_SIZE];
__complex128 c_l[MATRIX_SIZE][MATRIX_SIZE];
__complex128 c_u[MATRIX_SIZE][MATRIX_SIZE];
__complex128 c_y[MATRIX_SIZE];
__complex128 c_frb[MATRIX_SIZE];

__float128 f_f[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_r[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_fr[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_fra[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_l[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_u[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_y[MATRIX_SIZE];
__float128 f_frb[MATRIX_SIZE];

__float128 f_fa[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_fb[MATRIX_SIZE];
__float128 ap[MATRIX_SIZE][MATRIX_SIZE];
__float128 bp[MATRIX_SIZE];

// iteration
__complex128 c_w[MATRIX_SIZE];
__complex128 c_v[MATRIX_SIZE];
__complex128 c_z[MATRIX_SIZE];
__complex128 c_ax[MATRIX_SIZE];

__float128 f_w[MATRIX_SIZE];
__float128 f_v[MATRIX_SIZE];
__float128 f_z[MATRIX_SIZE];
__float128 f_ax[MATRIX_SIZE];


void solve_with_partial_pivot(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE]){
  int i,j;
  int p[MATRIX_SIZE]; // the last element is not used

  lu_partial_pivot_float128(a, f_l, f_u, p);

  for(i=0; i<MATRIX_SIZE; i++)
    bp[i] = b[i];
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      ap[i][j] = a[i][j];
  for(i=0; i<MATRIX_SIZE-1; i++){
    __float128 tmp = bp[i];
    bp[i] = bp[p[i]];
    bp[p[i]] = tmp;
  }
  for(i=0; i<MATRIX_SIZE-1; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      __float128 tmp = ap[i][j];
      ap[i][j] = ap[p[i]][j];
      ap[p[i]][j] = tmp;
    }
  }

  l_step_float128(f_l, bp, f_y);
  u_step_float128(f_u, f_y, x);

  for(i=0; i<MATRIX_SIZE; i++)
    xi[i] = x[i];

  iteration_float128(ap,f_l, f_u, bp, xi);
}

void iteration_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = 2;

  while(itera++ < iteramax){
    mat_vec_dot_float128(a,x,f_ax);
    vec_sub_float128(b,f_ax,f_w);
    l_step_float128(l,f_w,f_v);
    u_step_float128(u,f_v,f_z);
    vec_add_float128(x,f_z,x);
  }
}
void iteration_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 l[MATRIX_SIZE][MATRIX_SIZE], __complex128 u[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = 2;

  while(itera++ < iteramax){
    mat_vec_dot_complex128(a,x,c_ax);
    vec_sub_complex128(b,c_ax,c_w);
    l_step_complex128(l,c_w,c_v);
    u_step_complex128(u,c_v,c_z);
    vec_add_complex128(x,c_z,x);
  }
}

void solve_with_rdft_iteration(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE], __complex128 xi[MATRIX_SIZE]){
  int i;

  dft_matrix(c_f);
  r_matrix(c_r);
  mat_mul_complex128(c_f,c_r,c_fr);
  mat_mul_complex128(c_fr,a,c_fra);
  mat_vec_dot_complex128(c_fr,b,c_frb);

  lu_complex128(c_fra, c_l, c_u);

  l_step_complex128(c_l, c_frb, c_y);
  u_step_complex128(c_u, c_y, x);

  for(i=0; i<MATRIX_SIZE; i++)
    xi[i] = x[i];

  iteration_complex128(c_fra, c_l, c_u, c_frb, xi);
}
void solve_with_rdht_iteration(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE]){
  int i;

  dht_matrix(f_f);
  r_real_matrix(f_r);
  mat_mul_float128(f_f,f_r,f_fr);
  mat_mul_float128(f_fr,a,f_fra);
  mat_vec_dot_float128(f_fr,b,f_frb);

  lu_float128(f_fra, f_l, f_u);

  l_step_float128(f_l, f_frb, f_y);
  u_step_float128(f_u, f_y, x);

  for(i=0; i<MATRIX_SIZE; i++)
    xi[i] = x[i];

  iteration_float128(f_fra, f_l, f_u, f_frb, xi);
}

void solve_with_gauss_iteration(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE]){
  int i;

  gauss_matrix(f_f);
  mat_mul_float128(f_f,a,f_fa);
  mat_vec_dot_float128(f_f,b,f_fb);

  lu_float128(f_fa, f_l, f_u);

  l_step_float128(f_l, f_fb, f_y);
  u_step_float128(f_u, f_y, x);

  for(i=0; i<MATRIX_SIZE; i++)
    xi[i] = x[i];

  iteration_float128(f_fa, f_l, f_u, f_fb, xi);
}


