#ifndef _SCF_H_INCLUDED
#define _SCF_H_INCLUDED

#include <vector>

#include <Eigen/Dense>
#include <libint2.hpp>

#include "../lib/basis/BasisSet.h"
#include "../lib/chem/Molecule.h"
#include "../lib/wf/RhfWavefunction.h"
#include "../lib/wf/UhfWavefunction.h"

namespace minichem {

RhfWavefunction* rhf(Kernel* ker, BasisSet* bs, Molecule* mol);
UhfWavefunction* uhf(Kernel* ker, BasisSet* bs, Molecule* mol);

typedef std::vector<libint2::Shell> AtomCenteredBasis_t;

AtomCenteredBasis_t makeAtomCenteredSet(Molecule* mol, BasisSet* bs);
size_t basisDimension(const AtomCenteredBasis_t& shells);
size_t maxNumberPrimitives(const AtomCenteredBasis_t& shells);
int maxAngularMomentum(const AtomCenteredBasis_t& shells);
std::vector<size_t> mapShellBfn(const AtomCenteredBasis_t& shells);

#ifndef _MATRIX_TYPEDEF
#define _MATRIX_TYPEDEF
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
#endif

Matrix computeOneBodyInts(const AtomCenteredBasis_t& shells,
                          libint2::Operator obtype, Molecule* mol);

} // namespace minichem

#endif // _SCF_H_INCLUDED
