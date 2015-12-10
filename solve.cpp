#include "solve.h"
#include "lu.h"
#include "rdft.h"
#include "matlib.h"
#include <cstdio>

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

double d_f[MATRIX_SIZE][MATRIX_SIZE];
double d_r[MATRIX_SIZE][MATRIX_SIZE];
double d_fr[MATRIX_SIZE][MATRIX_SIZE];
double d_fra[MATRIX_SIZE][MATRIX_SIZE];
double d_l[MATRIX_SIZE][MATRIX_SIZE];
double d_u[MATRIX_SIZE][MATRIX_SIZE];
double d_y[MATRIX_SIZE];
double d_frb[MATRIX_SIZE];

std::complex<double> cd_f[MATRIX_SIZE][MATRIX_SIZE];
std::complex<double> cd_r[MATRIX_SIZE][MATRIX_SIZE];
std::complex<double> cd_fr[MATRIX_SIZE][MATRIX_SIZE];
std::complex<double> cd_fra[MATRIX_SIZE][MATRIX_SIZE];
std::complex<double> cd_l[MATRIX_SIZE][MATRIX_SIZE];
std::complex<double> cd_u[MATRIX_SIZE][MATRIX_SIZE];
std::complex<double> cd_y[MATRIX_SIZE];
std::complex<double> cd_frb[MATRIX_SIZE];

// partial pivot
double d_fa[MATRIX_SIZE][MATRIX_SIZE];
double d_fb[MATRIX_SIZE];
double d_ap[MATRIX_SIZE][MATRIX_SIZE];
double d_perm[MATRIX_SIZE][MATRIX_SIZE];
double d_bp[MATRIX_SIZE];

__float128 f_fa[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_fb[MATRIX_SIZE];
__float128 f_ap[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_perm[MATRIX_SIZE][MATRIX_SIZE];
__float128 f_bp[MATRIX_SIZE];

// iteration
__complex128 c_w[MATRIX_SIZE];
__complex128 c_v[MATRIX_SIZE];
__complex128 c_z[MATRIX_SIZE];
__complex128 c_ax[MATRIX_SIZE];

__float128 f_w[MATRIX_SIZE];
__float128 f_v[MATRIX_SIZE];
__float128 f_z[MATRIX_SIZE];
__float128 f_ax[MATRIX_SIZE];

double d_w[MATRIX_SIZE];
double d_v[MATRIX_SIZE];
double d_z[MATRIX_SIZE];
double d_ax[MATRIX_SIZE];

std::complex<double> cd_w[MATRIX_SIZE];
std::complex<double> cd_v[MATRIX_SIZE];
std::complex<double> cd_z[MATRIX_SIZE];
std::complex<double> cd_ax[MATRIX_SIZE];


double tmp1[MATRIX_SIZE][MATRIX_SIZE];
double tmp2[MATRIX_SIZE];
void solve_with_partial_pivot_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE], double xi[MATRIX_SIZE], double xia[MATRIX_SIZE]){
  int i,j;
  int p[MATRIX_SIZE]; // the last element is not used

  lu_partial_pivot_double(a, d_l, d_u, p);

  for(i=0; i<MATRIX_SIZE; i++)
    d_bp[i] = b[i];
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      d_ap[i][j] = a[i][j];
      d_perm[i][j] = 0.0;
    }
    d_perm[i][i] = 1.0;
  }
  for(i=0; i<MATRIX_SIZE-1; i++){
    double tmp = d_bp[i];
    d_bp[i] = d_bp[p[i]];
    d_bp[p[i]] = tmp;
  }
  for(i=0; i<MATRIX_SIZE-1; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      double tmp = d_ap[i][j];
      d_ap[i][j] = d_ap[p[i]][j];
      d_ap[p[i]][j] = tmp;
      tmp = d_perm[i][j];
      d_perm[i][j] = d_perm[p[i]][j];
      d_perm[p[i]][j] = tmp;
    }
  }

  l_step_double(d_l, d_bp, d_y);
  u_step_double(d_u, d_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_double(d_ap,d_l, d_u, d_bp, xi);
  iteration_double_another(d_perm, a,d_l, d_u, b, xia);
}
void solve_with_partial_pivot_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE], __float128 xia[MATRIX_SIZE]){
  int i,j;
  int p[MATRIX_SIZE]; // the last element is not used

  lu_partial_pivot_float128(a, f_l, f_u, p);

  for(i=0; i<MATRIX_SIZE; i++)
    f_bp[i] = b[i];
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      f_ap[i][j] = a[i][j];
      f_perm[i][j] = 0.0Q;
    }
    f_perm[i][i] = 1.0Q;
  }
  for(i=0; i<MATRIX_SIZE-1; i++){
    __float128 tmp = f_bp[i];
    f_bp[i] = f_bp[p[i]];
    f_bp[p[i]] = tmp;
  }
  for(i=0; i<MATRIX_SIZE-1; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      __float128 tmp = f_ap[i][j];
      f_ap[i][j] = f_ap[p[i]][j];
      f_ap[p[i]][j] = tmp;
      tmp = f_perm[i][j];
      f_perm[i][j] = f_perm[p[i]][j];
      f_perm[p[i]][j] = tmp;
    }
  }

  l_step_float128(f_l, f_bp, f_y);
  u_step_float128(f_u, f_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_float128(f_ap,f_l, f_u, f_bp, xi);
  iteration_float128_another(f_perm, a,f_l, f_u, b, xia);
}

void iteration_double(double a[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = 1;

  while(itera++ < iteramax){
    mat_vec_dot_double(a,x,d_ax);
    vec_sub_double(b,d_ax,d_w);
    l_step_double(l,d_w,d_v);
    u_step_double(u,d_v,d_z);
    vec_add_double(x,d_z,x);
  }
}
void iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = 1;

  while(itera++ < iteramax){
    mat_vec_dot_complex_double(a,x,cd_ax);
    vec_sub_complex_double(b,cd_ax,cd_w);
    l_step_complex_double(l,cd_w,cd_v);
    u_step_complex_double(u,cd_v,cd_z);
    vec_add_complex_double(x,cd_z,x);
  }
}
void iteration_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = 1;

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
  int iteramax = 1;

  while(itera++ < iteramax){
    mat_vec_dot_complex128(a,x,c_ax);
    vec_sub_complex128(b,c_ax,c_w);
    l_step_complex128(l,c_w,c_v);
    u_step_complex128(u,c_v,c_z);
    vec_add_complex128(x,c_z,x);
  }
}
void iteration_double_another(double f[MATRIX_SIZE][MATRIX_SIZE], double a[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = 1;

  while(itera++ < iteramax){
    mat_vec_dot_double(a,x,d_ax);
    vec_sub_double(b,d_ax,d_w);
    mat_vec_dot_double(f,d_w,d_ax);
    l_step_double(l,d_ax,d_v);
    u_step_double(u,d_v,d_z);
    vec_add_double(x,d_z,x);
  }
}
void iteration_complex_double_another(std::complex<double> f[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = 1;

  while(itera++ < iteramax){
    mat_vec_dot_complex_double(a,x,cd_ax);
    vec_sub_complex_double(b,cd_ax,cd_w);
    mat_vec_dot_complex_double(f,cd_w,cd_ax);
    l_step_complex_double(l,cd_ax,cd_v);
    u_step_complex_double(u,cd_v,cd_z);
    vec_add_complex_double(x,cd_z,x);
  }
}
void iteration_float128_another(__float128 f[MATRIX_SIZE][MATRIX_SIZE], __float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 l[MATRIX_SIZE][MATRIX_SIZE], __float128 u[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = 1;

  while(itera++ < iteramax){
    mat_vec_dot_float128(a,x,f_ax);
    vec_sub_float128(b,f_ax,f_w);
    mat_vec_dot_float128(f,f_w,f_ax);
    l_step_float128(l,f_ax,f_v);
    u_step_float128(u,f_v,f_z);
    vec_add_float128(x,f_z,x);
  }
}
void iteration_complex128_another(__complex128 f[MATRIX_SIZE][MATRIX_SIZE], __complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 l[MATRIX_SIZE][MATRIX_SIZE], __complex128 u[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = 1;

  while(itera++ < iteramax){
    mat_vec_dot_complex128(a,x,c_ax);
    vec_sub_complex128(b,c_ax,c_w);
    mat_vec_dot_complex128(f,c_w,c_ax);
    l_step_complex128(l,c_ax,c_v);
    u_step_complex128(u,c_v,c_z);
    vec_add_complex128(x,c_z,x);
  }
}

void solve_with_rdft_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE]){
  int i;

  dft_matrix_complex_double(cd_f);
  r_matrix_complex_double(cd_r);
  mat_mul_complex_double(cd_f,cd_r,cd_fr);
  mat_mul_complex_double(cd_fr,a,cd_fra);
  mat_vec_dot_complex_double(cd_fr,b,cd_frb);

  lu_complex_double(cd_fra, cd_l, cd_u);

  l_step_complex_double(cd_l, cd_frb, cd_y);
  u_step_complex_double(cd_u, cd_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_complex_double(cd_fra, cd_l, cd_u, cd_frb, xi);
  iteration_complex_double_another(cd_fr, a, cd_l, cd_u, b, xia);
}
void solve_with_rdft_iteration_complex128(__complex128 a[MATRIX_SIZE][MATRIX_SIZE], __complex128 b[MATRIX_SIZE], __complex128 x[MATRIX_SIZE], __complex128 xi[MATRIX_SIZE], __complex128 xia[MATRIX_SIZE]){
  int i;

  dft_matrix_complex_128(c_f);
  r_matrix_complex_128(c_r);
  mat_mul_complex128(c_f,c_r,c_fr);
  mat_mul_complex128(c_fr,a,c_fra);
  mat_vec_dot_complex128(c_fr,b,c_frb);

  lu_complex128(c_fra, c_l, c_u);

  l_step_complex128(c_l, c_frb, c_y);
  u_step_complex128(c_u, c_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_complex128(c_fra, c_l, c_u, c_frb, xi);
  iteration_complex128_another(c_fr, a, c_l, c_u, b, xia);
}
void solve_with_rdht_iteration_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE], __float128 xia[MATRIX_SIZE]){
  int i;

  dht_matrix(f_f);
  r_real_matrix(f_r);
  mat_mul_float128(f_f,f_r,f_fr);
  mat_mul_float128(f_fr,a,f_fra);
  mat_vec_dot_float128(f_fr,b,f_frb);

  lu_float128(f_fra, f_l, f_u);

  l_step_float128(f_l, f_frb, f_y);
  u_step_float128(f_u, f_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_float128(f_fra, f_l, f_u, f_frb, xi);
  iteration_float128_another(f_fr, a, f_l, f_u, b, xia);
}

void solve_with_gauss_iteration_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE], double xi[MATRIX_SIZE], double xia[MATRIX_SIZE]){
  int i;

  gauss_matrix_double(d_f);
  mat_mul_double(d_f,a,d_fa);
  mat_vec_dot_double(d_f,b,d_fb);

  lu_double(d_fa, d_l, d_u);

  l_step_double(d_l, d_fb, d_y);
  u_step_double(d_u, d_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_double(d_fa, d_l, d_u, d_fb, xi);
  iteration_double_another(d_f, a, d_l, d_u, b, xia);
}
void solve_with_gauss_iteration_float128(__float128 a[MATRIX_SIZE][MATRIX_SIZE], __float128 b[MATRIX_SIZE], __float128 x[MATRIX_SIZE], __float128 xi[MATRIX_SIZE], __float128 xia[MATRIX_SIZE]){
  int i;

  gauss_matrix_float128(f_f);
  mat_mul_float128(f_f,a,f_fa);
  mat_vec_dot_float128(f_f,b,f_fb);

  lu_float128(f_fa, f_l, f_u);

  l_step_float128(f_l, f_fb, f_y);
  u_step_float128(f_u, f_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_float128(f_fa, f_l, f_u, f_fb, xi);
  iteration_float128_another(f_f, a, f_l, f_u, b, xia);
}


