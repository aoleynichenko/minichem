/* Restricted Hartree-Fock Method
 *
 * Note: this code is based on the "hartree-fock" example from the Libint
 * library by Edward Valeev and coworkers. For original code, see
 * https://github.com/evaleev/libint/tree/master/tests/hartree-fock
 */

#include <algorithm>
#include <array>
#include <sstream>
#include <vector>

#include "../Kernel.h"
#include "../lib/basis/BasisSet.h"
#include "../lib/chem/Molecule.h"
#include "../lib/wf/RhfWavefunction.h"

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>

#include "scf.h"

namespace minichem {

using libint2::Engine;
using libint2::Operator;
using libint2::Shell;
using std::ostringstream;
using std::vector;

OutputStream* out = nullptr;

AtomCenteredBasis_t makeAtomCenteredSet(Molecule* mol, BasisSet* bs);
size_t basisDimension(const AtomCenteredBasis_t& shells);
size_t maxNumberPrimitives(const AtomCenteredBasis_t& shells);
int maxAngularMomentum(const AtomCenteredBasis_t& shells);
std::vector<size_t> mapShellBfn(const AtomCenteredBasis_t& shells);

Matrix computeOneBodyInts(const AtomCenteredBasis_t& shells,
                          Operator obtype, Molecule* mol);

RhfWavefunction* rhf(Kernel* ker, BasisSet* bs, Molecule* mol)
{
  BasisSet basis = bs->filter(mol);
  out = ker->getOutput();

  out->printf("               ********************************\n");
  out->printf("               *      Hartree-Fock Method     *\n");
  out->printf("               ********************************\n");
  out->println();
  out->printf("                     Basis Set Information\n");
  out->printf("                    -----------------------\n");
  out->printf("%s", basis.toString().c_str());

  // basis set: create and print information
  auto shells = makeAtomCenteredSet(mol, bs);
  size_t nao = basisDimension(shells);

  out->printf("Number of shells   : %d\n", shells.size());
  out->printf("Dimension          : %d\n", nao);
  out->printf("Integration engine : %s\n", "Libint");

  // initializes the Libint integrals library ... now ready to compute
  libint2::initialize();

  Matrix S = computeOneBodyInts(shells, Operator::overlap, mol); // overlap
  Matrix T = computeOneBodyInts(shells, Operator::kinetic, mol); // kinetic
  Matrix V = computeOneBodyInts(shells, Operator::nuclear, mol); // nuclear-attraction
  Matrix H = T + V; // core Hamiltonian
  // T and V no longer needed, free up the memory
  T.resize(0,0);
  V.resize(0,0);
  mainlog->log("One-electron integrals done");

  // AO basis set orthogonalization

  // initial guess

  // main loop
  libint2::finalize(); // done with libint

  return nullptr;
}

AtomCenteredBasis_t makeAtomCenteredSet(Molecule* mol, BasisSet* bs)
{
  AtomCenteredBasis_t bfns;
  vector<Atom>* atoms = mol->getAtoms();

  for (auto a = atoms->begin(); a != atoms->end(); a++) {
    BasisSet::BasisTemplate_t* blocks = bs->getBlocks(a->charge);
    for (auto b = blocks->begin(); b != blocks->end(); b++) // loop over L-blocks
      for (auto c = b->contr_.begin(); c != b->contr_.end(); c++) { // loop over contracted f-ns
        // one-by-one! Generally contracted basis sets are yet unavailable in Libint
        vector<Shell::Contraction> libint_contr; // libint-style contracted f-ns
        libint_contr.push_back({b->l_, !b->cart_, *c});
        // !!! here: implicit contruction of libint2::Shell object from
        // vector<double>, vector<Contraction> and std::array<double, 3>
        bfns.push_back({b->alpha_, libint_contr, {a->x, a->y, a->z}});
      }
  }
  return bfns;
}

Matrix computeOneBodyInts(const AtomCenteredBasis_t& shells,
                          libint2::Operator obtype,
                          Molecule* mol)
{
  const auto n = basisDimension(shells);
  auto max_l = maxAngularMomentum(shells);  // is const?
  auto max_nprim = maxNumberPrimitives(shells);
  auto atoms = mol->getAtoms();
  Matrix result(n, n);

  // construct the 1-e integrals engine
  // Engine can be native for minichem (not Libint's)
  Engine engine(obtype, max_nprim, max_l, 0);

  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  // of course, point nuclei are not suitable for relativistic treatment!
  if (obtype == Operator::nuclear) {
    vector<std::pair<double,std::array<double,3>>> q;
    for (auto a = atoms->begin(); a != atoms->end(); a++)
      q.push_back({static_cast<double>(a->charge), {a->x, a->y, a->z}});
    engine.set_params(q);
  }

  auto offsets = mapShellBfn(shells);

  // buf[0] points to the target shell set
  // after every call to engine.compute()
  const auto& buf = engine.results();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
  // we COMPUTE rectangular BLOCK of the 'results' matrix!
  for(auto s1 = 0; s1 != shells.size(); s1++) {
    auto bf1 = offsets[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();
    for(auto s2 = 0; s2 <= s1; s2++) {
      auto bf2 = offsets[s2];
      auto n2 = shells[s2].size();
      // compute shell pair; return is the pointer to the buffer
      engine.compute(shells[s1], shells[s2]);
      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
      // type buf_mat == Eigen::Matrix
      Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
    }
  }

  return result;
}

// utilities for working with AtomCenteredBasis_t
// based on code from Libint example hartree-fock
size_t basisDimension(const AtomCenteredBasis_t& shells)
{
  size_t n = 0;
  for (const auto& shell: shells)
    n += shell.size();
  return n;
}

size_t maxNumberPrimitives(const AtomCenteredBasis_t& shells)
{
  size_t n = 0;
  for (auto shell: shells)
    n = std::max(shell.nprim(), n);
  return n;
}

int maxAngularMomentum(const AtomCenteredBasis_t& shells)
{
  int l = 0;
  for (auto shell: shells)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}

// returns offsets of shells in 1-dim array [ao-bfn1, ..., ao-bfnN]
vector<size_t> mapShellBfn(const AtomCenteredBasis_t& shells)
{
  vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto shell: shells) {
    result.push_back(n);
    n += shell.size();
  }
  return result;
}

} // namespace minichem
