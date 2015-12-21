#include "givens.h"
#include "gen.h"
#include <cmath>
#include <vector>
#include <random>

extern int for_perm[MATRIX_SIZE];
extern int init_for_perm_flg;
extern void init_for_perm();

std::complex<double> givens_tmp_vector[MATRIX_SIZE];
std::complex<double> givens_tmp_matrix[MATRIX_SIZE][MATRIX_SIZE];

// ----- create/delete givens matrix basic ----- //
givens_matrix *create_givens_matrix(int i, int j, double theta){
  givens_matrix *gm = new givens_matrix;

  if(i > j){
    int tmp = i;
    i = j;
    j = tmp;
  }

  gm->i = i;
  gm->j = j;
  gm->c = std::complex<double>(cos(theta),0);
  gm->s = std::complex<double>(sin(theta),0);

  return gm;
}
void delete_givens_matrix(givens_matrix *gm){
  if(gm) delete gm;
}

givens_matrix_list *create_givens_matrix_list(){
  givens_matrix_list *gml = new givens_matrix_list;
  gml->gm = NULL;
  gml->next = NULL;
  return gml;
}
void delete_givens_matrix_list(givens_matrix_list *gml){
  if(gml){
    if(gml->next != NULL){
      delete_givens_matrix_list(gml->next);
    }
    delete_givens_matrix(gml->gm);
    delete gml;
  }
}

// ----- create/delete givens matrix ----- //
givens_matrix_list *r_givens_matrix_double(){
  givens_matrix_list *gml_r, *gml_cur;
  std::random_device rd;
	std::mt19937 mt(rd());

  if(! init_for_perm_flg){
    init_for_perm();
  }
  std::vector<int> perm(for_perm,for_perm+MATRIX_SIZE);

  gml_r = NULL;
  gml_cur = gml_r;
  while(perm.size() > 1){
    int p,q;
    int idx = mt() % (perm.size());
    p = perm[idx];
    perm.erase(perm.begin()+idx);
    idx = mt() % (perm.size());
    q = perm[idx];
    perm.erase(perm.begin()+idx);

    double rval = uniform()*M_PI;
    if(rval < 0.0)
      rval = -rval;
    givens_matrix_list *gml_tmp = create_givens_matrix_list();
    givens_matrix *gm_tmp = create_givens_matrix(p,q,rval);
    gml_tmp->gm   = gm_tmp;
    gml_tmp->next = NULL;

    if(gml_r == NULL){
      gml_r = gml_tmp;
      gml_cur = gml_r;
    }else{
      gml_cur->next = gml_tmp;
      gml_cur = gml_cur->next;
    }
  }
  return gml_r;
}
std::complex<double> qq[MATRIX_SIZE][MATRIX_SIZE];
std::complex<double> pp[MATRIX_SIZE][MATRIX_SIZE];
void ins_r_givens_matrix_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE]){
  std::random_device rd;
	std::mt19937 mt(rd());

  if(! init_for_perm_flg){
    init_for_perm();
  }
  std::vector<int> perm(for_perm,for_perm+MATRIX_SIZE);

  for(int i=0; i<MATRIX_SIZE; i++){
    for(int j=0; j<MATRIX_SIZE; j++){
      qq[i][j] = std::complex<double>(0.0,0.0);
    }
    qq[i][i] = std::complex<double>(1.0,0.0);
  }
  while(perm.size() > 1){
    int p,q;
    int idx = mt() % (perm.size());
    p = perm[idx];
    perm.erase(perm.begin()+idx);
    idx = mt() % (perm.size());
    q = perm[idx];
    perm.erase(perm.begin()+idx);

    double rval = uniform()*M_PI;
    if(rval < 0.0)
      rval = -rval;
    for(int i=0; i<MATRIX_SIZE; i++){
      for(int j=0; j<MATRIX_SIZE; j++){
        pp[i][j] = std::complex<double>(0.0,0.0);
      }
      pp[i][i] = std::complex<double>(1.0,0.0);
    }
    pp[p][p] = std::complex<double>(cos(rval),0);
    pp[q][p] = std::complex<double>(-sin(rval),0);
    pp[p][q] = std::complex<double>(sin(rval),0);
    pp[q][q] = std::complex<double>(cos(rval),0);
    mat_mul_complex_double(qq,pp,a);
    for(int i=0; i<MATRIX_SIZE; i++){
      for(int j=0; j<MATRIX_SIZE; j++){
        qq[i][j] = a[i][j];
      }
    }
  }
}


// ----- convert givens matrix ----- //
void convert_givens_matrix_complex_double(givens_matrix *gm, std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE]){
    int i,j;
    for(i=0; i<MATRIX_SIZE; i++){
      for(j=0; j<MATRIX_SIZE; j++){
        a[i][j] = std::complex<double>(0.0,0.0);
      }
      a[i][i]= std::complex<double>(1.0,0.0);
    }
    a[gm->i][gm->i] = gm->c;
    a[gm->i][gm->j] = -(gm->s);
    a[gm->j][gm->i] = gm->s;
    a[gm->j][gm->j] = gm->c;
}

// ----- givens matrix-vector multiplication ----- //
void mat_vec_dot_givens_complex_double_(givens_matrix *gm, std::complex<double> a[MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE]){
  int i;
  for(i=0; i<MATRIX_SIZE; i++){
    if(i == gm->i){
      b[i] = (gm->c)*a[gm->i] + (-(gm->s))*a[gm->j];
    }else if(i == gm->j){
      b[i] = (gm->s)*a[gm->i] + (gm->c)*a[gm->j];
    }else{
      b[i] = a[i];
    }
  }
}

void mat_vec_dot_givens_complex_double(givens_matrix_list *gml, std::complex<double> a[MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE]){
  int i;
  givens_matrix_list *gml_cur = gml;
  for(i=0; i<MATRIX_SIZE; i++)
    givens_tmp_vector[i] = a[i];
  while(gml_cur != NULL){
    mat_vec_dot_givens_complex_double_(gml_cur->gm, givens_tmp_vector, b);
    gml_cur = gml_cur->next;
    for(i=0; i<MATRIX_SIZE; i++)
      givens_tmp_vector[i] = b[i];
  }
}

// ----- givens matrix-matrix multiplication ----- //
//void mat_mul_givens_left_complex_double_(givens_matrix *gm, std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE]){
//  int i,j;
//  for(i=0; i<MATRIX_SIZE; i++){
//    for(j=0; j<MATRIX_SIZE; j++){
//      if(i != gm->i || i != gm->j){
//        b[i][j] = a[i][j];
//      }else if(i == gm->i){
//        b[i][j] = (gm->c)*a[i][j] + (-(gm->s))*a[j][j];
//      }else{ // i == gm->j
//        b[i][j] = (gm->s)*a[i][j] + (gm->c)*a[j][j];
//      }
//    }
//  }
//}
void mat_mul_givens_right_complex_double_(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], givens_matrix *gm, std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      if(j == gm->i){
        b[i][gm->i] = a[i][gm->i]*(gm->c) + a[i][gm->j]*(gm->s);
      }else if(j == gm->j){
        b[i][gm->j] = a[i][gm->i]*(-(gm->s)) + a[i][gm->j]*(gm->c);
      }else{
        b[i][j] = a[i][j];
      }
    }
  }
}
void mat_mul_givens_left_complex_double(givens_matrix_list *gml, std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE]){
  givens_matrix_list *gml_cur = gml;
  while(gml_cur != NULL){
    //mat_mul_givens_left_complex_double_(gml_cur->gm, a, b);
    gml_cur = gml_cur->next;
  }
}
void mat_mul_givens_right_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], givens_matrix_list *gml, std::complex<double> b[MATRIX_SIZE][MATRIX_SIZE]){
  int i,j;
  givens_matrix_list *gml_cur = gml;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      givens_tmp_matrix[i][j] = a[i][j];
  while(gml_cur != NULL){
    mat_mul_givens_right_complex_double_(givens_tmp_matrix, gml_cur->gm, b);
    gml_cur = gml_cur->next;
    for(i=0; i<MATRIX_SIZE; i++)
      for(j=0; j<MATRIX_SIZE; j++)
        givens_tmp_matrix[i][j] = b[i][j];
  }
}
