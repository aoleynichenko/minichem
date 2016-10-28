#include "../Kernel.h"
#include "../lib/basis/BasisSet.h"
#include "../lib/chem/Molecule.h"
#include "../lib/wf/RhfWavefunction.h"

#include "scf.h"

namespace minichem {

RhfWavefunction* rhf(Kernel* ker, BasisSet* bs, Molecule* mol)
{
  BasisSet basis = bs->filter(mol);
  OutputStream* out = ker->getOutput();

  out->printf("               ********************************\n");
  out->printf("               *      Hartree-Fock Method     *\n");
  out->printf("               ********************************\n");
  out->println();
  out->printf("%s\n", basis.toString().c_str());

  return nullptr;
}

} // namespace minichem
