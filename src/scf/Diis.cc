#include "Diis.h"
#include <iostream>
using namespace std;

#include <Eigen/Dense>

namespace minichem {

Diis::Diis(int diisbas, const Matrix& S)
  : bas_m(diisbas), S_m(S)
{
}

// returns max error in ErrMatrix
double Diis::storeFock(const Matrix& F, const Matrix& D)
{
  Matrix E = F*D*S_m - S_m*D*F;
  listFE_m.push_back({F, E});
  if (listFE_m.size() > bas_m)
    listFE_m.pop_front();
  return E.maxCoeff();
}

// Frobenius matrix product
double mdot(const Matrix& A, const Matrix& B)
{
  auto nrows = A.rows();
  auto ncols = A.cols();
  double prod = 0.0;
  for (auto i = 0; i < nrows; i++)
    for (auto j = 0; j < ncols; j++)
      prod += A(i,j) * B(i,j);
  return prod;
}

Matrix Diis::extrapolate()
{
  auto diislen = listFE_m.size();

  Matrix B = Matrix::Zero(diislen+1,diislen+1);
  for (auto i = 0; i < diislen; i++) {
		for (auto j = i; j < diislen; j++) {
			B(i,j) = mdot(listFE_m[i].second, listFE_m[j].second);  // B_ij = e_i*e_j
      B(j,i) = B(i,j);
    }
    B(i,diislen) = -1.0;
	}
	// -1 -1 -1 ... 0
	for (auto i = 0; i < diislen; i++)
		B(diislen,i) = -1.0;
	B(diislen,diislen) = 0.0;

  // A*x = B  --->  B*c = (0,0,0...-1) = r = right (is a column!)
  Eigen::VectorXd right(diislen+1);
	for (auto i = 0; i < diislen; i++)
		right[i] = 0.0;
	right[diislen] = -1;

  // solve system of linear equations
  Eigen::VectorXd c = B.colPivHouseholderQr().solve(right);

  // extrapole Fock matrix
  // TODO: no error check here!
  auto fdim = listFE_m[0].first.rows();
  Matrix F = Matrix::Zero(fdim, fdim);
  for (auto i = 0; i < diislen; i++)
    F += c[i] * listFE_m[i].first;

  return F;
}

Diis::~Diis()
{
}

} //namespace minichem
