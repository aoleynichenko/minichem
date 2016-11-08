/* Unrestricted Hartree-Fock Method
 */

#include <algorithm>
#include <array>
#include <sstream>
#include <vector>

#include "../Kernel.h"
#include "../lib/basis/BasisSet.h"
#include "../lib/chem/Molecule.h"
#include "../lib/wf/UhfWavefunction.h"

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

static OutputStream* out = nullptr;

Matrix computeTwoBodyPart_uhf(const std::vector<libint2::Shell>& shells,
                          const Matrix& D);

UhfWavefunction* uhf(Kernel* ker, BasisSet* bs, Molecule* mol)
{
  // timing
  struct {
    double time_1e;
    double time_guess;
    double time_fock;
    double time_diag;
    double time_dens;
    double time_diis;
  } uhf_timing = {0, 0, 0, 0, 0, 0};

  BasisSet basis = bs->filter(mol);
  out = ker->getOutput();

  out->printf("               ******************************************\n");
  out->printf("               *    Unrestricted Hartree-Fock Method    *\n");
  out->printf("               ******************************************\n");
  out->println();
  out->printf("                     Basis Set Information\n");
  out->printf("                    -----------------------\n");
  out->printf("%s", basis.toString().c_str());

  // basis set: create and print information
  auto shells = makeAtomCenteredSet(mol, bs);
  size_t nao = basisDimension(shells);
  auto na = mol->nalpha();
  auto nb = mol->nbeta();

  out->printf("Nalpha             : %d\n", na);
  out->printf("Nbeta              : %d\n", nb);
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
  uhf_timing.time_1e += std::chrono::duration_cast<std::chrono::milliseconds>(t1 - start).count();
  mainlog->log("One-electron integrals done");

  // initial guess
  Matrix Da, Db;
  t1 = std::chrono::system_clock::now();
  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
  auto Ca = gen_eig_solver.eigenvectors();
  auto Cb = Ca;
  auto Ca_occ = Ca.leftCols(na);
  Da = Ca_occ * Ca_occ.transpose();
  auto Cb_occ = Cb.leftCols(nb);
  Db = Cb_occ * Cb_occ.transpose();
  auto t2 = std::chrono::system_clock::now();
  uhf_timing.time_guess += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

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
    auto Da_last = Da, Db_last = Db;

    // build a new Fock matrix
    auto t3 = std::chrono::system_clock::now();
    auto Fa = H, Fb = H;
    //F += computeTwoBodyPart_uhf(shells, D);
    auto t4 = std::chrono::system_clock::now();
    uhf_timing.time_fock += std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();

    // DIIS
    /*if (iter > 2) {
      diis.storeFock(F, D);
      if (iter > 3)
        F = diis.extrapolate(); // obtain new F from DIIS extrapolation
    }
    auto t_diis = std::chrono::system_clock::now();
    rhf_timing.time_diis += std::chrono::duration_cast<std::chrono::milliseconds>(t_diis - t4).count();
    */

    // solve F C = e S C
    auto t5 = std::chrono::system_clock::now();
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver_a(Fa, S);
    auto Ca = gen_eig_solver.eigenvectors();
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver_b(Fb, S);
    auto Cb = gen_eig_solver.eigenvectors();
    auto t6 = std::chrono::system_clock::now();
    uhf_timing.time_diag += std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();

    // compute density, D = C(occ) . C(occ)T
    auto t7 = std::chrono::system_clock::now();
    auto Ca_occ = Ca.leftCols(na);
    Da = Ca_occ * Ca_occ.transpose();
    auto Cb_occ = Cb.leftCols(nb);
    Db = Cb_occ * Cb_occ.transpose();
    auto t8 = std::chrono::system_clock::now();
    uhf_timing.time_dens += std::chrono::duration_cast<std::chrono::milliseconds>(t8 - t7).count();

    // compute HF energy
    ehf = 0.0;
    for (auto i = 0; i < nao; i++)
      for (auto j = 0; j < nao; j++) {
        ehf += 0.5 * ((Da(i,j) + Db(i,j))*H(j,i) + Da(j,i)*Fa(i,j) + Db(j,i)*Fb(i,j));
      }

    // compute difference with last iteration
    ediff = ehf - ehf_last;
    rmsd = (Da+Db - Da_last - Db_last).norm();  // computed using total density matrix

    auto t = std::chrono::system_clock::now();
    curr_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t -start).count();
    out->printf(" %02d %20.12f %20.12f %20.12f %20.12f %13.3f\n", iter, ehf, ehf + enuc, ediff, rmsd, curr_ms/1000);
  } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

  out->printf("----------------------------------------------------------------------------------------------------------\n");
  out->printf(" Converged!\n\n");
  out->printf("           Total UHF energy = %20.12f\n", ehf + enuc);
  out->printf("          Electronic energy = %20.12f\n", ehf);
//  out->printf("        One-electron energy = %20.12f\n", e1e);
//  out->printf("        Two-electron energy = %20.12f\n", e2e);
  out->printf("   Nuclear repulsion energy = %20.12f\n", enuc);

  // print timing
  out->println();
  out->printf("      Time for:           sec\n");
  out->printf("---------------------------------\n");
  out->printf("  One-electron ints    %9.3f\n", uhf_timing.time_1e/1000);
  out->printf("  Initial guess        %9.3f\n", uhf_timing.time_guess/1000);
  out->printf("  Fock matrix          %9.3f\n", uhf_timing.time_fock/1000);
  out->printf("  Density matrix       %9.3f\n", uhf_timing.time_dens/1000);
  out->printf("  Diagonalization      %9.3f\n", uhf_timing.time_diag/1000);
  out->printf("  DIIS extrapolation   %9.3f\n", uhf_timing.time_diis/1000);
  out->printf("  Time per iteration   %9.3f\n", curr_ms/iter/1000);
  out->printf("---------------------------------\n\n");

  libint2::finalize(); // done with libint

  return nullptr;
}

Matrix computeTwoBodyPart_uhf(const std::vector<libint2::Shell>& shells,
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

} // namespace minichem
