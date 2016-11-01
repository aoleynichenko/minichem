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
Matrix computeTwoBodyPart_simple(const AtomCenteredBasis_t& shells,
                          const Matrix& D);

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
  auto nocc = mol->nelec()/2;

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

  // initial guess
  Matrix D;
  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
  auto C = gen_eig_solver.eigenvectors();
  auto C_occ = C.leftCols(nocc);
  D = C_occ * C_occ.transpose();

  // main loop
  const auto maxiter = 100;
  const auto conv = 1e-12;
  auto iter = 0;
  auto rmsd = 0.0;
  auto ediff = 0.0;
  auto ehf = 0.0, e1e = 0.0, e2e = 0.0;
  auto enuc = mol->enuc();

  out->printf(" Iter        E(elec)              E(tot)               Delta(E)             RMS(D)    \n");
  out->printf("--------------------------------------------------------------------------------------\n");
  do {
    iter++;

    // Save a copy of the energy and the density
    auto ehf_last = ehf;
    auto D_last = D;

    // build a new Fock matrix
    auto F = H;
    F += computeTwoBodyPart_simple(shells, D);

    // solve F C = e S C
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
    auto C = gen_eig_solver.eigenvectors();

    // compute density, D = C(occ) . C(occ)T
    auto C_occ = C.leftCols(nocc);
    D = C_occ * C_occ.transpose();

    // compute HF energy
    e1e = e2e = 0.0;
    for (auto i = 0; i < nao; i++)
      for (auto j = 0; j < nao; j++) {
        e1e += D(i,j) * H(i,j);
        e2e += D(i,j) * F(i,j);
	}
	ehf = e1e + e2e;

    // compute difference with last iteration
    ediff = ehf - ehf_last;
    rmsd = (D - D_last).norm();

    out->printf(" %02d %20.12f %20.12f %20.12f %20.12f\n", iter, ehf, ehf + enuc, ediff, rmsd);
  } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));
  
  out->printf("--------------------------------------------------------------------------------------\n");
  out->printf(" Converged!\n\n");
  out->printf("           Total RHF energy = %20.12f\n", ehf + enuc);
  out->printf("          Electronic energy = %20.12f\n", ehf);
  out->printf("        One-electron energy = %20.12f\n", e1e);
  out->printf("        Two-electron energy = %20.12f\n", e2e);
  out->printf("   Nuclear repulsion energy = %20.12f\n", enuc);

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
        libint_contr.push_back({b->l_, b->cart_, *c});
        // !!! here: implicit contruction of libint2::Shell object from
        // vector<double>, vector<Contraction> and std::array<double, 3>
        bfns.push_back({b->alpha_, libint_contr, {{a->x, a->y, a->z}}});
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

// the simplest (and the slowest!) implementation - without permutational
// symmetry of 2-e integrals
Matrix computeTwoBodyPart_simple(const AtomCenteredBasis_t& shells,
                                 const Matrix& D)
{
  const auto n = basisDimension(shells);
  auto max_l = maxAngularMomentum(shells);  // is const?
  auto max_nprim = maxNumberPrimitives(shells);
  Matrix G = Matrix::Zero(n,n);

  // construct the ERI engine
  Engine engine(Operator::coulomb, max_nprim, max_l, 0);
  auto offs = mapShellBfn(shells);
  const auto& buf = engine.results();

  // loop over shell pairs of the Fock matrix, {s1,s2}
  // Fock matrix is symmetric, but skipping it here for simplicity
  for(auto s1 = 0; s1 != shells.size(); s1++) {
    auto bf1_first = offs[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();
    for(auto s2 = 0; s2 != shells.size(); s2++) {
      auto bf2_first = offs[s2];
      auto n2 = shells[s2].size();

      // loop over shell pairs of the density matrix, {s3,s4}
      // again symmetry is not used for simplicity
      for(auto s3 = 0; s3 != shells.size(); s3++) {
        auto bf3_first = offs[s3];
        auto n3 = shells[s3].size();
        for(auto s4 = 0; s4 != shells.size(); s4++) {
          auto bf4_first = offs[s4];
          auto n4 = shells[s4].size();

          // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
          // Compute shell quartet
          // J = <PQ|RS>    K = <PQ|SR>
          // J = (PR|QS)    K = (PS|QR)
          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* eris_J = buf[0];
          if (eris_J == nullptr)
            continue; // if all integrals screened out, skip to next quartet

          // hence some manual labor here:
          // 1) loop over every integral in the shell set (= nested loops over basis functions in each shell)
          // and 2) add contribution from each integral
          for(auto f1 = 0, idx = 0; f1 != n1; f1++) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2 = 0; f2 != n2; f2++) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3 = 0; f3 != n3; f3++) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4 = 0; f4 != n4; ++f4, ++idx) {  // f4++ ??
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) += D(bf3,bf4) * 2.0 * eris_J[idx];
                }
              }
            }
          }

          // exchange contribution to the Fock matrix is from {s1,s3,s2,s4} integrals
          engine.compute(shells[s1], shells[s3], shells[s2], shells[s4]);
          const auto* eris_K = buf[0];

          for(auto f1 = 0, idx = 0; f1 != n1; f1++) {
            const auto bf1 = f1 + bf1_first;
            for(auto f3 = 0; f3 != n3; f3++) {
              const auto bf3 = f3 + bf3_first;
              for(auto f2 = 0; f2 != n2; f2++) {
                const auto bf2 = f2 + bf2_first;
                for(auto f4 = 0; f4 != n4; ++f4, ++idx) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) -= D(bf3,bf4) * eris_K[idx];
                }
              }
            }
          }

        }
      }
    }
  }

  return G;
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
