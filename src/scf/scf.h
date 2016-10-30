#ifndef _SCF_H_INCLUDED
#define _SCF_H_INCLUDED

#include <vector>

#include <Eigen/Dense>
#include <libint2.hpp>

#include "../lib/basis/BasisSet.h"
#include "../lib/chem/Molecule.h"
#include "../lib/wf/RhfWavefunction.h"

namespace minichem {

RhfWavefunction* rhf(Kernel* ker, BasisSet* bs, Molecule* mol);

typedef std::vector<libint2::Shell> AtomCenteredBasis_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

} // namespace minichem

#endif // _SCF_H_INCLUDED
