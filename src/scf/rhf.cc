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
#include "Diis.h"

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
Matrix computeTwoBodyPart(const std::vector<libint2::Shell>& shells,
                          const Matrix& D);

RhfWavefunction* rhf(Kernel* ker, BasisSet* bs, Molecule* mol)
{
  // timing
  struct {
    double time_1e;
    double time_guess;
    double time_fock;
    double time_diag;
    double time_dens;
    double time_diis;
  } rhf_timing = {0, 0, 0, 0, 0};

  BasisSet basis = bs->filter(mol);
  out = ker->getOutput();

  out->printf("               ******************************************\n");
  out->printf("               *     Restricted Hartree-Fock Method     *\n");
  out->printf("               ******************************************\n");
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
  out->printf("Initial guess      : %s\n", "Hcore");
  out->println();

  // initializes the Libint integrals library ... now ready to compute
  libint2::initialize();

  // one-electron matrices
  auto start = std::chrono::system_clock::now();
  Matrix S = computeOneBodyInts(shells, Operator::overlap, mol); // overlap
  Matrix T = computeOneBodyInts(shells, Operator::kinetic, mol); // kinetic
  Matrix V = computeOneBodyInts(shells, Operator::nuclear, mol); // nuclear-attraction
  Matrix H = T + V; // core Hamiltonian
  // T and V no longer needed, free up the memory
  T.resize(0,0);
  V.resize(0,0);
  auto t1 = std::chrono::system_clock::now();
  rhf_timing.time_1e += std::chrono::duration_cast<std::chrono::milliseconds>(t1 - start).count();
  mainlog->log("One-electron integrals done");

  // initial guess
  Matrix D;
  t1 = std::chrono::system_clock::now();
  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
  auto C = gen_eig_solver.eigenvectors();
  auto C_occ = C.leftCols(nocc);
  D = C_occ * C_occ.transpose();
  auto t2 = std::chrono::system_clock::now();
  rhf_timing.time_guess += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  // main loop
  const auto maxiter = 150;
  const auto conv = 1e-8;
  auto iter = 0;
  auto rmsd = 0.0;
  auto ediff = 0.0;
  auto ehf = 0.0, e1e = 0.0, e2e = 0.0;
  auto enuc = mol->enuc();
  double curr_ms = 0.0;
  Diis diis(6, S);

  out->printf(" Iter        E(elec)              E(tot)               Delta(E)              RMS(D)          Time, sec\n");
  out->printf("----------------------------------------------------------------------------------------------------------\n");
  do {
    iter++;

    // Save a copy of the energy and the density
    auto ehf_last = ehf;
    auto D_last = D;

    // build a new Fock matrix
    auto t3 = std::chrono::system_clock::now();
    auto F = H;
    F += computeTwoBodyPart(shells, D);
    auto t4 = std::chrono::system_clock::now();
    rhf_timing.time_fock += std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();

    // DIIS
    if (iter > 2) {
      diis.storeFock(F, D);
      if (iter > 3)
        F = diis.extrapolate(); // obtain new F from DIIS extrapolation
    }
    auto t_diis = std::chrono::system_clock::now();
    rhf_timing.time_diis += std::chrono::duration_cast<std::chrono::milliseconds>(t_diis - t4).count();

    // solve F C = e S C
    auto t5 = std::chrono::system_clock::now();
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
    auto C = gen_eig_solver.eigenvectors();
    auto t6 = std::chrono::system_clock::now();
    rhf_timing.time_diag += std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();

    // compute density, D = C(occ) . C(occ)T
    auto t7 = std::chrono::system_clock::now();
    auto C_occ = C.leftCols(nocc);
    D = C_occ * C_occ.transpose();
    auto t8 = std::chrono::system_clock::now();
    rhf_timing.time_dens += std::chrono::duration_cast<std::chrono::milliseconds>(t8 - t7).count();

    // compute HF energy
    e1e = e2e = 0.0;
    for (auto i = 0; i < nao; i++)
      for (auto j = 0; j < nao; j++) {
        e1e += D(i,j) * H(i,j) * 2;         // One-electron energy
        e2e += D(i,j) * (F(i,j) - H(i,j));  // Two-electron energy (in fact, * G(i,j))
      }
    ehf = e1e + e2e;

    // compute difference with last iteration
    ediff = ehf - ehf_last;
    rmsd = (D - D_last).norm();

    auto t = std::chrono::system_clock::now();
    curr_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t -start).count();
    out->printf(" %02d %20.12f %20.12f %20.12f %20.12f %13.3f\n", iter, ehf, ehf + enuc, ediff, rmsd, curr_ms/1000);
  } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

  out->printf("----------------------------------------------------------------------------------------------------------\n");
  out->printf(" Converged!\n\n");
  out->printf("           Total RHF energy = %20.12f\n", ehf + enuc);
  out->printf("          Electronic energy = %20.12f\n", ehf);
  out->printf("        One-electron energy = %20.12f\n", e1e);
  out->printf("        Two-electron energy = %20.12f\n", e2e);
  out->printf("   Nuclear repulsion energy = %20.12f\n", enuc);

  // print timing
  out->println();
  out->printf("      Time for:           sec\n");
  out->printf("---------------------------------\n");
  out->printf("  One-electron ints    %9.3f\n", rhf_timing.time_1e/1000);
  out->printf("  Initial guess        %9.3f\n", rhf_timing.time_guess/1000);
  out->printf("  Fock matrix          %9.3f\n", rhf_timing.time_fock/1000);
  out->printf("  Density matrix       %9.3f\n", rhf_timing.time_dens/1000);
  out->printf("  Diagonalization      %9.3f\n", rhf_timing.time_diag/1000);
  out->printf("  DIIS extrapolation   %9.3f\n", rhf_timing.time_diis/1000);
  out->printf("  Time per iteration   %9.3f\n", curr_ms/iter/1000);
  out->printf("---------------------------------\n\n");

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

Matrix computeTwoBodyPart(const std::vector<libint2::Shell>& shells,
                          const Matrix& D)
{
  const auto n = basisDimension(shells);
  auto max_l = maxAngularMomentum(shells);  // is const?
  auto max_nprim = maxNumberPrimitives(shells);
  Matrix G = Matrix::Zero(n,n);

  // construct the 2-electron repulsion integrals engine
  Engine engine(Operator::coulomb, max_nprim, max_l, 0);

  auto shell2bf = mapShellBfn(shells);

  const auto& buf = engine.results();

  // The problem with the simple Fock builder is that permutational symmetries of the Fock,
  // density, and two-electron integrals are not taken into account to reduce the cost.
  // To make the simple Fock builder efficient we must rearrange our computation.
  // The most expensive step in Fock matrix construction is the evaluation of 2-e integrals;
  // hence we must minimize the number of computed integrals by taking advantage of their permutational
  // symmetry. Due to the multiplicative and Hermitian nature of the Coulomb kernel (and realness
  // of the Gaussians) the permutational symmetry of the 2-e ints is given by the following relations:
  //
  // (12|34) = (21|34) = (12|43) = (21|43) = (34|12) = (43|12) = (34|21) = (43|21)
  //
  // (here we use chemists' notation for the integrals, i.e in (ab|cd) a and b correspond to
  // electron 1, and c and d -- to electron 2).
  //
  // It is easy to verify that the following set of nested loops produces a permutationally-unique
  // set of integrals:
  // foreach a = 0 .. n-1
  //   foreach b = 0 .. a
  //     foreach c = 0 .. a
  //       foreach d = 0 .. (a == c ? b : c)
  //         compute (ab|cd)
  //
  // The only complication is that we must compute integrals over shells. But it's not that complicated ...
  //
  // The real trick is figuring out to which matrix elements of the Fock matrix each permutationally-unique
  // (ab|cd) contributes. STOP READING and try to figure it out yourself. (to check your answer see below)

  // loop over permutationally-unique set of shells
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();   // number of basis functions in this shell

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      for(auto s3=0; s3<=s1; ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        const auto s4_max = (s1 == s3) ? s2 : s3;
        for(auto s4=0; s4<=s4_max; ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
          auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
          auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue; // if all integrals screened out, skip to next quartet

          // ANSWER
          // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
          //    F(a,b) += (ab|cd) * D(c,d)
          //    F(c,d) += (ab|cd) * D(a,b)
          //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
          //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
          //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
          //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
          // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
          //    i.e. the number of the integrals/sets equivalent to it
          // 3) the end result must be symmetrized
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;

                  const auto value = buf_1234[f1234];

                  const auto value_scal_by_deg = value * s1234_deg;

                  G(bf1,bf2) += D(bf3,bf4) * value_scal_by_deg;
                  G(bf3,bf4) += D(bf1,bf2) * value_scal_by_deg;
                  G(bf1,bf3) -= 0.25 * D(bf2,bf4) * value_scal_by_deg;
                  G(bf2,bf4) -= 0.25 * D(bf1,bf3) * value_scal_by_deg;
                  G(bf1,bf4) -= 0.25 * D(bf2,bf3) * value_scal_by_deg;
                  G(bf2,bf3) -= 0.25 * D(bf1,bf4) * value_scal_by_deg;
                }
              }
            }
          }

        }
      }
    }
  }

  // symmetrize the result and return
  Matrix Gt = G.transpose();
  return 0.5 * (G + Gt);
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
