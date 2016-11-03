#ifndef _DIIS_H_INCLUDED
#define _DIIS_H_INCLUDED

#include <deque>
#include <utility>

#include <Eigen/Core>

namespace minichem {

#ifndef _MATRIX_TYPEDEF
#define _MATRIX_TYPEDEF
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
#endif

// RHF-specific DIIS
// Maybe it is a good idea to make this class a template?
// Diis<SCF>, Diis<CC> ... m?
class Diis {
public:
  Diis(int diisbas, const Matrix& S); // S is needed for E = FDS-SDF
  double storeFock(const Matrix& F, const Matrix& D); // returns max error in ErrMatrix
  Matrix extrapolate();
  ~Diis();
private:
  int bas_m;
  const Matrix& S_m;
  std::deque<std::pair<Matrix,Matrix>> listFE_m;
};

} // namespace minichem

#endif // _DIIS_H_INCLUDED
