#ifndef _SCF_H_INCLUDED
#define _SCF_H_INCLUDED

#include "../lib/basis/BasisSet.h"
#include "../lib/chem/Molecule.h"
#include "../lib/wf/RhfWavefunction.h"

namespace minichem {

RhfWavefunction* rhf(Kernel* ker, BasisSet* bs, Molecule* mol);

} // namespace minichem

#endif // _SCF_H_INCLUDED
