#include "eigen.h"

#include <Eigen/Core>
#include <Eigen/SVD>

using namespace Eigen;

double condition_number(double a[MATRIX_SIZE][MATRIX_SIZE]){
  MatrixXd mat =  Map<Matrix<double, Dynamic, Dynamic> >(&(a[0][0]), MATRIX_SIZE, MATRIX_SIZE);
  JacobiSVD<MatrixXd> svd(mat);

  auto sval = svd.singularValues();

  return sval[0]/sval[sval.size()-1];
}
