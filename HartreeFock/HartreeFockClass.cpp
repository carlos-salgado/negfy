/*
 * HartreeFockClass.cpp
 *
 *  Created on: 18 mar. 2017
 *      Author: carlos
 */





#include "HartreeFockClass.h"
/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <thread>
#include <atomic>
#include <iterator>
#include <unordered_map>
#include <mutex>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string.h> // to use stricmp.

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
//#include <Eigen/IOFormat> // IS A CLASS.

// have BTAS library?
#ifdef LIBINT2_HAVE_BTAS
//# include <btas/btas.h>
#endif // LIBINT2_HAVE_BTAS


// Libint Gaussian integrals library
#include <libint2/diis.h>
#include <libint2/util/intpart_iter.h>
#include <libint2/chemistry/sto3g_atomic_density.h>
#include <libint2.hpp>

//#include <libint2/shell.h>
//#include <libint2/atom.h>
//#include <libint2/basis.h>


#if defined(_OPENMP)
//# include <omp.h>
#endif




//****************************************************************************************************************************************
//************************************** DFT PART FROM ERKALE ****************************************************************************
//****************************************************************************************************************************************


/* COMENTADA DOBLEMENTE, DEJAR SOLO LO DE ABAJO.
#include "global.h"
#include <algorithm>
#include <armadillo>
#include <cstdio>
#include <cfloat>
#include "basis.h"
#include "broyden.h"
#include "elements.h"
#include "dftfuncs.h"
#include "dftgrid.h"
#include "diis.h"
#include "guess.h"
#include "linalg.h"
#include "mathf.h"
#include "properties.h"
//#include "scf.h" // TIENE LOS STRUCT COMO SOL, EN... y las declaraciones de void Fock_UDFT(...);
#include "stringutil.h"
#include "timer.h"
#include "trrh.h"
extern "C" {
#include <gsl/gsl_poly.h>
}
*/


#include "settings.h"
#include "checkpoint.h"
#include "erkalebasis.h"
#include "dftgrid.h"
#include "dftfuncs.h"
#include "scf.h"
#include "trrh.h"
#include <algorithm>
#include <armadillo>
#include <cstdio>
#include <cfloat>


#include <vector>
//using ErkaleSCF::SCF;
//using SCF;
//using DFTGrid;


//****************************************************************************************************************************************
//************************************ END DFT PART FROM ERKALE **************************************************************************
//****************************************************************************************************************************************








using libint2::Shell;
using libint2::Atom;
using libint2::BasisSet;
using libint2::Operator;
using libint2::BraKet;

using namespace std;

namespace libint2 {
int nthreads; 
}


namespace HartreeFock {
/*
extern "C" {
    void matrix_multiply(double *A, int *rows_A, int *cols_A, double *B, int *rows_B, int *cols_B, double *C, int *rows_C, int *cols_C);
};
*/

HartreeFockClass::HartreeFockClass() {
	// TODO Auto-generated constructor stub

}

HartreeFockClass::~HartreeFockClass() {
	// TODO Auto-generated destructor stub
}

//Fill a vector with random numbers in the range [lower, upper]
void HartreeFockClass::rnd_fill(vector<double> &V, double lower, double upper) {

    //use system clock to create a random seed
    size_t seed (clock());

    //use the default random engine and an uniform distribution
    default_random_engine eng(seed);
    uniform_real_distribution<double> distr(lower, upper);

    for( auto &elem : V){
        elem = distr(eng);
    }
}

//Print matrix V(rows, cols) storage in column-major format
void HartreeFockClass::print_matrix(vector<double> const &V, int rows, int cols) {

    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
            cout << V[j * rows + i] << " ";
        }
        cout << std::endl;
    }
    cout << std::endl;
}

/*
void mat_mult(double *A, int *rows_A, int *cols_A, double *B, int *rows_B, int *cols_B, double *C, int *rows_C, int *cols_C){
	//matrix_multiply(*A, *rows_A, *cols_A, *B, *rows_B, *cols_B, *C, *rows_C, *cols_C);
}
*/
//void mat_mult(double *A, int *rows_A, int *cols_A, double *B, int *rows_B, int *cols_B, double *C, int *rows_C, int *cols_C){
//void HartreeFockClass::mat_mult(std::vector<double, std::allocator<double>> A, int rows_A, int cols_A, std::vector<double, std::allocator<double>> B, int rows_B, int cols_B, std::vector<double, std::allocator<double>> C, int rows_C, int cols_C){
//matrix_multiply(*A, *rows_A, *cols_A, *B, *rows_B, *cols_B, *C, *rows_C, *cols_C);
//void HartreeFockClass::mat_mult(vector<double> const &A, int rows_A, int cols_A,vector<double> const &B, int rows_B, int cols_B,vector<double> const &C, int rows_C, int cols_C){


/* 2018-04-30 FROM FORTRAN MODULE MATRIXMULTIPLY, ELIMINATED.
void HartreeFockClass::mat_mult(vector<double> &A, int rows_A, int cols_A,vector<double> &B, int rows_B, int cols_B,vector<double> &C, int rows_C, int cols_C){
	HartreeFock::matrix_multiply(&A[0], &rows_A, &cols_A, &B[0], &rows_B, &cols_B, &C[0], &rows_C, &cols_C);
}
*/






//int main(int argc, char *argv[]) {
void HartreeFockClass::run(char *geom, char *basis) {
  using std::cout;
  using std::cerr;
  using std::endl;

  try {

    /*** =========================== ***/
    /*** initialize molecule         ***/
    /*** =========================== ***/

    // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
//    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
//    const auto basisname = (argc > 2) ? argv[2] : "aug-cc-pVDZ";
//    const auto filename = (argc > 1) ? geom : "h2o.xyz";
//    const auto basisname = (argc > 2) ? basis : "aug-cc-pVDZ";
    const auto filename = "h2o.xyz";
    //const auto basisname = "aug-cc-pVDZ";
    const auto basisname = "STO-3G";

    bool do_density_fitting = false;
#ifdef HAVE_DENSITY_FITTING
    do_density_fitting = (argc > 3);
    const auto dfbasisname = do_density_fitting ? argv[3] : "";
#endif
    std::vector<Atom> atoms = read_geometry(filename);

    // set up thread pool
    {
      using libint2::nthreads;
      auto nthreads_cstr = getenv("LIBINT_NUM_THREADS");
      nthreads = 1;
      if (nthreads_cstr && strcmp(nthreads_cstr, "")) {
        std::istringstream iss(nthreads_cstr);
        iss >> nthreads;
        if (nthreads > 1 << 16 || nthreads <= 0) nthreads = 1;
      }
#if defined(_OPENMP)
      omp_set_num_threads(nthreads);
#endif
      cout << "Will scale over " << nthreads
#if defined(_OPENMP)
                << " OpenMP"
#else
                << " C++11"
#endif
                << " threads" << std::endl;
    }

    // count the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i) nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;
    cout << "# of electrons = " << nelectron << endl;

    // compute the nuclear repulsion energy
    auto enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij * xij + yij * yij + zij * zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    cout << "Nuclear repulsion energy = " << std::setprecision(15) << enuc
         << endl;

    libint2::Shell::do_enforce_unit_normalization(false);

    cout << "Atomic Cartesian coordinates (a.u.):" << endl;
    for (const auto& a : atoms)
      cout << a.atomic_number << " " << a.x << " " << a.y << " " << a.z
                << std::endl;

    BasisSet obs(basisname, atoms);
    cout << "orbital basis set rank = " << obs.nbf() << endl;

#ifdef HAVE_DENSITY_FITTING
    BasisSet dfbs;
    if (do_density_fitting) {
      dfbs = BasisSet(dfbasisname, atoms);
      cout << "density-fitting basis set rank = " << dfbs.nbf() << endl;
    }
#endif  // HAVE_DENSITY_FITTING

    /*** =========================== ***/
    /*** compute 1-e integrals       ***/
    /*** =========================== ***/

    // initializes the Libint integrals library ... now ready to compute
    libint2::initialize();

    // compute OBS non-negligible shell-pair list
    {
      obs_shellpair_list = compute_shellpair_list(obs);
      size_t nsp = 0;
      for (auto& sp : obs_shellpair_list) {
        nsp += sp.second.size();
      }
      cout << "# of {all,non-negligible} shell-pairs = {"
                << obs.size() * (obs.size() + 1) / 2 << "," << nsp << "}"
                << std::endl;
    }

    // compute one-body integrals
    auto S = compute_1body_ints<Operator::overlap>(obs)[0];
    auto T = compute_1body_ints<Operator::kinetic>(obs)[0];
    auto V = compute_1body_ints<Operator::nuclear>(obs, atoms)[0];
    Matrix H = T + V;
    T.resize(0, 0);
    V.resize(0, 0);

    // compute orthogonalizer X such that X.transpose().S.X = I
    Matrix X, Xinv;
    double XtX_condition_number;  // condition number of "re-conditioned"
                                  // overlap obtained as Xinv.transpose().Xinv
    // one should think of columns of Xinv as the conditioned basis
    // Re: name ... cond # (Xinv.transpose().Xinv) = cond # (X.transpose() .
    // X)
    // by default assume can manage to compute with condition number of S <=
    // 1/eps
    // this is probably too optimistic, but in well-behaved cases even 10^11 is
    // OK
    double S_condition_number_threshold =
        1.0 / std::numeric_limits<double>::epsilon();
    std::tie(X, Xinv, XtX_condition_number) =
        conditioning_orthogonalizer(S, S_condition_number_threshold);

    Matrix D;
    Matrix C_occ;
    Matrix evals;
    {  // use SOAD as the guess density
      const auto tstart = std::chrono::high_resolution_clock::now();

      auto D_minbs = compute_soad(atoms);  // compute guess in minimal basis
      //BasisSet minbs("STO-3G", atoms);
      //"aug-cc-pVDZ"
      BasisSet minbs(basisname, atoms);

      if (minbs == obs)
        D = D_minbs;
      else {  // if basis != minimal basis, map non-representable SOAD guess
              // into the AO basis
              // by diagonalizing a Fock matrix
        cout << "projecting SOAD into AO basis ... ";
        auto F = H;
        F += compute_2body_fock_general(
            obs, D_minbs, minbs, true /* SOAD_D_is_shelldiagonal */,
            std::numeric_limits<double>::epsilon()  // this is cheap, no reason
                                                    // to be cheaper
            );

        // solve F C = e S C by (conditioned) transformation to F' C' = e C',
        // where
        // F' = X.transpose().F.X; the original C is obtained as C = X.C'
        Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(X.transpose() * F * X);
        auto C = X * eig_solver.eigenvectors();

        // compute density, D = C(occ).C(occ)T
        C_occ = C.leftCols(ndocc);
        D = C_occ * C_occ.transpose();

        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;
        cout << "done (" << time_elapsed.count() << " s)" << std::endl;
      }
    }

    // pre-compute data for Schwarz bounds
    auto K = compute_schwarz_ints<>(obs);

// prepare for density fitting
#ifdef HAVE_DENSITY_FITTING
    std::unique_ptr<DFFockEngine> dffockengine(
        do_density_fitting ? new DFFockEngine(obs, dfbs) : nullptr);
#endif  // HAVE_DENSITY_FITTING

    /*** =========================== ***/
    /***          SCF loop           ***/
    /*** =========================== ***/

    const auto maxiter = 100;
    const auto conv = 1e-12;
    auto iter = 0;
    auto rms_error = 1.0;
    auto ediff_rel = 0.0;
    auto ehf = 0.0;
    auto n2 = D.cols() * D.rows();
    libint2::DIIS<Matrix> diis(2);  // start DIIS on second iteration

    // prepare for incremental Fock build ...
    Matrix D_diff = D;
    Matrix F = H;
    bool reset_incremental_fock_formation = false;
    bool incremental_Fbuild_started = false;
    double start_incremental_F_threshold = 1e-5;
    double next_reset_threshold = 0.0;
    size_t last_reset_iteration = 0;
    // ... unless doing DF, then use MO coefficients, hence not "incremental"
    if (do_density_fitting) start_incremental_F_threshold = 0.0;

    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;

      // Last iteration's energy and density
      auto ehf_last = ehf;
      Matrix D_last = D;

      if (not incremental_Fbuild_started &&
          rms_error < start_incremental_F_threshold) {
        incremental_Fbuild_started = true;
        reset_incremental_fock_formation = false;
        last_reset_iteration = iter - 1;
        next_reset_threshold = rms_error / 1e1;
        cout << "== started incremental fock build" << std::endl;
      }
      if (reset_incremental_fock_formation || not incremental_Fbuild_started) {
        F = H;
        D_diff = D;
      }
      if (reset_incremental_fock_formation && incremental_Fbuild_started) {
        reset_incremental_fock_formation = false;
        last_reset_iteration = iter;
        next_reset_threshold = rms_error / 1e1;
        cout << "== reset incremental fock build" << std::endl;
      }

      // build a new Fock matrix
      if (not do_density_fitting) {
        // totally empirical precision variation, involves the condition number
        const auto precision_F = std::min(
            std::min(1e-3 / XtX_condition_number, 1e-7),
            std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon()));
        F += compute_2body_fock(obs, D_diff, precision_F, K);
      }
#if HAVE_DENSITY_FITTING
      else {  // do DF
        F = H + dffockengine->compute_2body_fock_dfC(C_occ);
      }
#else
      else {
        assert(false);
      }  // do_density_fitting is true but HAVE_DENSITY_FITTING is not defined!
         // should not happen
#endif  // HAVE_DENSITY_FITTING

      // compute HF energy with the non-extrapolated Fock matrix
      ehf = D.cwiseProduct(H + F).sum();
      ediff_rel = std::abs((ehf - ehf_last) / ehf);

      // compute SCF error
      Matrix FD_comm = F * D * S - S * D * F;
      rms_error = FD_comm.norm() / n2;
      if (rms_error < next_reset_threshold || iter - last_reset_iteration >= 8)
        reset_incremental_fock_formation = true;

      // DIIS extrapolate F
      Matrix F_diis = F;  // extrapolated F cannot be used in incremental Fock
                          // build; only used to produce the density
                          // make a copy of the unextrapolated matrix
      diis.extrapolate(F_diis, FD_comm);

      // solve F C = e S C by (conditioned) transformation to F' C' = e C',
      // where
      // F' = X.transpose().F.X; the original C is obtained as C = X.C'
      Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(X.transpose() * F_diis *
                                                       X);
      evals = eig_solver.eigenvalues();
      auto C = X * eig_solver.eigenvectors();

      // compute density, D = C(occ).C(occ)T
      C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();
      D_diff = D - D_last;

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;

      if (iter == 1)
        cout << "\nIter         E(HF)                 D(E)/E         "
                     "RMS([F,D])/nn       Time(s)\n";
      printf(" %02d %20.12f %20.12e %20.12e %10.5lf\n", iter, ehf + enuc,
             ediff_rel, rms_error, time_elapsed.count());

    } while (((ediff_rel > conv) || (rms_error > conv)) && (iter < maxiter));

    printf("** Hartree-Fock energy = %20.12f\n", ehf + enuc);

    auto Mu = compute_1body_ints<Operator::emultipole2>(obs);
    std::array<double, 3> mu;
    std::array<double, 6> qu;
    for (int xyz = 0; xyz != 3; ++xyz)
      mu[xyz] = -2 *
                D.cwiseProduct(Mu[xyz + 1])
                    .sum();  // 2 = alpha + beta, -1 = electron charge
    for (int k = 0; k != 6; ++k)
      qu[k] = -2 *
              D.cwiseProduct(Mu[k + 4])
                  .sum();  // 2 = alpha + beta, -1 = electron charge
    cout << "** edipole = ";
    std::copy(mu.begin(), mu.end(),
              std::ostream_iterator<double>(cout, " "));
    cout << std::endl;
    cout << "** equadrupole = ";
    std::copy(qu.begin(), qu.end(),
              std::ostream_iterator<double>(cout, " "));
    cout << std::endl;

    {  // compute force
#if LIBINT2_DERIV_ONEBODY_ORDER
      // compute 1-e forces
      Matrix F1 = Matrix::Zero(atoms.size(), 3);
      Matrix F_Pulay = Matrix::Zero(atoms.size(), 3);
      //////////
      // one-body contributions to the forces
      //////////
      auto T1 = compute_1body_ints_deriv<Operator::kinetic>(1, obs, atoms);
      auto V1 = compute_1body_ints_deriv<Operator::nuclear>(1, obs, atoms);
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          auto force = 2 * (T1[i] + V1[i]).cwiseProduct(D).sum();
          F1(atom, xyz) += force;
        }
      }

      //////////
      // Pulay force
      //////////
      // orbital energy density
      DiagonalMatrix evals_occ(evals.topRows(ndocc));
      Matrix W = C_occ * evals_occ * C_occ.transpose();
      auto S1 = compute_1body_ints_deriv<Operator::overlap>(1, obs, atoms);
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          auto force = 2 * S1[i].cwiseProduct(W).sum();
          F_Pulay(atom, xyz) -= force;
        }
      }

      cout << "** 1-body forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << F1(atom, xyz) << " ";
      cout << std::endl;
      cout << "** Pulay forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz)
          cout << F_Pulay(atom, xyz) << " ";
      cout << std::endl;
#endif  // LIBINT2_DERIV_ONEBODY_ORDER

#if LIBINT2_DERIV_ERI_ORDER
      // compute 2-e forces
      Matrix F2 = Matrix::Zero(atoms.size(), 3);

      //////////
      // two-body contributions to the forces
      //////////
      auto G1 = compute_2body_fock_deriv<1>(obs, atoms, D);
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
          auto force = G1[i].cwiseProduct(D).sum();
          F2(atom, xyz) += force;
        }
      }

      cout << "** 2-body forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << F2(atom, xyz) << " ";
      cout << std::endl;
#endif

// if support 1-e and 2-e derivatives compute nuclear repulsion force and the
// total force
#if LIBINT2_DERIV_ONEBODY_ORDER && LIBINT2_DERIV_ERI_ORDER
      // compute nuclear repulsion forces
      Matrix FN = Matrix::Zero(atoms.size(), 3);
      //////////
      // nuclear repulsion contribution to the forces
      //////////
      for (auto a1 = 1; a1 != atoms.size(); ++a1) {
        const auto& atom1 = atoms[a1];
        for (auto a2 = 0; a2 < a1; ++a2) {
          const auto& atom2 = atoms[a2];

          auto x12 = atom1.x - atom2.x;
          auto y12 = atom1.y - atom2.y;
          auto z12 = atom1.z - atom2.z;
          auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
          auto r12 = sqrt(r12_2);
          auto r12_3 = r12 * r12_2;

          auto z1z2_over_r12_3 =
              atom1.atomic_number * atom2.atomic_number / r12_3;

          auto fx = -x12 * z1z2_over_r12_3;
          auto fy = -y12 * z1z2_over_r12_3;
          auto fz = -z12 * z1z2_over_r12_3;
          FN(a1, 0) += fx;
          FN(a1, 1) += fy;
          FN(a1, 2) += fz;
          FN(a2, 0) -= fx;
          FN(a2, 1) -= fy;
          FN(a2, 2) -= fz;
        }
      }

      cout << "** nuclear repulsion forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << FN(atom, xyz) << " ";
      cout << std::endl;

      auto F = F1 + F_Pulay + F2 + FN;
      cout << "** Hartree-Fock forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << F(atom, xyz) << " ";
      cout << std::endl;
#endif
    }

    {  // compute hessian
      const auto ncoords = 3 * atoms.size();
      // # of elems in upper triangle
      const auto nelem =  ncoords * (ncoords+1) / 2;
#if LIBINT2_DERIV_ONEBODY_ORDER > 1
      // compute 1-e hessian
      Matrix H1 = Matrix::Zero(ncoords, ncoords);
      Matrix H_Pulay = Matrix::Zero(ncoords, ncoords);
      //////////
      // one-body contributions to the hessian
      //////////
      auto T2 = compute_1body_ints_deriv<Operator::kinetic>(2, obs, atoms);
      auto V2 = compute_1body_ints_deriv<Operator::nuclear>(2, obs, atoms);
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          auto hess = 2 * (T2[i] + V2[i]).cwiseProduct(D).sum();
          H1(row, col) += hess;
        }
      }

      //////////
      // Pulay hessian
      //////////
      // orbital energy density
      DiagonalMatrix evals_occ(evals.topRows(ndocc));
      Matrix W = C_occ * evals_occ * C_occ.transpose();
      auto S2 = compute_1body_ints_deriv<Operator::overlap>(2, obs, atoms);
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          auto hess = 2 * S2[i].cwiseProduct(W).sum();
          H_Pulay(row, col) -= hess;
        }
      }

      cout << "** 1-body hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << H1(row, col) << " ";
        }
      }
      cout << std::endl;

      cout << "** Pulay hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << H_Pulay(row, col) << " ";
        }
      }
      cout << std::endl;
#endif  // LIBINT2_DERIV_ONEBODY_ORDER > 1

#if LIBINT2_DERIV_ERI_ORDER > 1
      // compute 2-e forces
      Matrix H2 = Matrix::Zero(ncoords, ncoords);

      //////////
      // two-body contributions to the forces
      //////////
      auto G2 = compute_2body_fock_deriv<2>(obs, atoms, D);
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
          auto hess = G2[i].cwiseProduct(D).sum();
          H2(row, col) += hess;
        }
      }

      cout << "** 2-body hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << H2(row, col) << " ";
        }
      }
      cout << std::endl;
#endif

// if support 1-e and 2-e 2nd derivatives compute nuclear repulsion hessian and
// the total hessian
#if LIBINT2_DERIV_ONEBODY_ORDER > 1 && LIBINT2_DERIV_ERI_ORDER > 1
      // compute nuclear repulsion hessian
      // NB only the upper triangle is computed!!!
      Matrix HN = Matrix::Zero(ncoords, ncoords);

      //////////
      // nuclear repulsion contribution to the hessian
      //////////
      for (auto a1 = 1; a1 != atoms.size(); ++a1) {
        const auto& atom1 = atoms[a1];
        for (auto a2 = 0; a2 < a1; ++a2) {
          const auto& atom2 = atoms[a2];

          auto x12 = atom1.x - atom2.x;
          auto y12 = atom1.y - atom2.y;
          auto z12 = atom1.z - atom2.z;
          auto x12_2 = x12 * x12;
          auto y12_2 = y12 * y12;
          auto z12_2 = z12 * z12;
          auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
          auto r12 = sqrt(r12_2);
          auto r12_5 = r12 * r12_2 * r12_2;

          auto z1z2_over_r12_5 =
              atom1.atomic_number * atom2.atomic_number / r12_5;

          HN(3*a1 + 0, 3*a1 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
          HN(3*a1 + 1, 3*a1 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
          HN(3*a1 + 2, 3*a1 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
          HN(3*a1 + 0, 3*a1 + 1) += z1z2_over_r12_5 * (3*x12*y12);
          HN(3*a1 + 0, 3*a1 + 2) += z1z2_over_r12_5 * (3*x12*z12);
          HN(3*a1 + 1, 3*a1 + 2) += z1z2_over_r12_5 * (3*y12*z12);

          HN(3*a2 + 0, 3*a2 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
          HN(3*a2 + 1, 3*a2 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
          HN(3*a2 + 2, 3*a2 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
          HN(3*a2 + 0, 3*a2 + 1) += z1z2_over_r12_5 * (3*x12*y12);
          HN(3*a2 + 0, 3*a2 + 2) += z1z2_over_r12_5 * (3*x12*z12);
          HN(3*a2 + 1, 3*a2 + 2) += z1z2_over_r12_5 * (3*y12*z12);

          HN(3*a2 + 0, 3*a1 + 0) -= z1z2_over_r12_5 * (3*x12_2 - r12_2);
          HN(3*a2 + 1, 3*a1 + 1) -= z1z2_over_r12_5 * (3*y12_2 - r12_2);
          HN(3*a2 + 2, 3*a1 + 2) -= z1z2_over_r12_5 * (3*z12_2 - r12_2);
          HN(3*a2 + 1, 3*a1 + 0) -= z1z2_over_r12_5 * (3*y12*x12);
          HN(3*a2 + 2, 3*a1 + 0) -= z1z2_over_r12_5 * (3*z12*x12);
          HN(3*a2 + 2, 3*a1 + 1) -= z1z2_over_r12_5 * (3*z12*y12);
          HN(3*a2 + 0, 3*a1 + 1) -= z1z2_over_r12_5 * (3*x12*y12);
          HN(3*a2 + 0, 3*a1 + 2) -= z1z2_over_r12_5 * (3*x12*z12);
          HN(3*a2 + 1, 3*a1 + 2) -= z1z2_over_r12_5 * (3*y12*z12);
        }
      }

      cout << "** nuclear repulsion hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << HN(row, col) << " ";
        }
      }
      cout << std::endl;

      auto H = H1 + H_Pulay + H2 + HN;
      cout << "** Hartree-Fock hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << H(row, col) << " ";
        }
      }
      cout << std::endl;
#endif
    }

    libint2::finalize();  // done with libint

  }  // end of try block; if any exceptions occurred, report them and exit
     // cleanly

  catch (const char* ex) {
    cerr << "caught exception: " << ex << endl;
    //return 1;
  } catch (std::string& ex) {
    cerr << "caught exception: " << ex << endl;
    //return 1;
  } catch (std::exception& ex) {
    cerr << ex.what() << endl;
    //return 1;
  } catch (...) {
    cerr << "caught unknown exception\n";
    //return 1;
  }

  //return 0;
}
//-------------------------------------------------------------------------------------------
//--------------------------- RUN ANT -------------------------------------------------------
//-------------------------------------------------------------------------------------------
void HartreeFockClass::runant(const char *basename, const char *basisname, const char *fullname, const bool debugflag, const bool erkaleflag) {
  using std::cout;
  using std::cerr;
  using std::endl;

  if (debugflag){
	  printf("-------------------------------------------------------------------------\n");
	  printf("--------------------------- DEBUGFLAG IS TRUE ---------------------------\n");
	  printf("-------------------------------------------------------------------------\n");
  }

  try {

    /*** =========================== ***/
    /*** initialize molecule         ***/
    /*** =========================== ***/

    // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
//    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
//    const auto basisname = (argc > 2) ? argv[2] : "aug-cc-pVDZ";
//    const auto filename = (argc > 1) ? geom : "h2o.xyz";
//    const auto basisname = (argc > 2) ? basis : "aug-cc-pVDZ";

	// THE LINE BELOW WORKS
	//const auto basename = "h2";
	char filename[strlen(basename)+strlen(basisname)+4];
	strcpy(filename,basename);
	//strcat(filename,"_");
	//strcat(filename,basisname);
	strcat(filename,".xyz");

	cout << "filename = "<<filename<<"\n";

	//const auto filename = "li2_STO-3G.xyz";
    //const auto basisname = "aug-cc-pVDZ";
	//const auto basisname = "augmentation-cc-pVDZ";

	// THE LINE BELOW WORKS
	//const auto basisname = "STO-3G";

	//bool do_density_fitting = false;
	//bool do_density_fitting = true; // SET TRUE TO PERFORM DFT WITH ERKALE.
	bool do_density_fitting = erkaleflag;

//#ifdef HAVE_DENSITY_FITTING // COMMENTED BY C.SALGADO ON 2019-04-30 TO AVOID ERRORS WITH DFT-HF
//    do_density_fitting = (argc > 3);
//    const auto dfbasisname = do_density_fitting ? argv[3] : "";
//#endif

    std::vector<Atom> atoms = read_geometry(filename);

    // set up thread pool
    {
      using libint2::nthreads;
      auto nthreads_cstr = getenv("LIBINT_NUM_THREADS");
      nthreads = 1;
      if (nthreads_cstr && strcmp(nthreads_cstr, "")) {
        std::istringstream iss(nthreads_cstr);
        iss >> nthreads;
        if (nthreads > 1 << 16 || nthreads <= 0) nthreads = 1;
      }
#if defined(_OPENMP)
      omp_set_num_threads(nthreads);
#endif
      cout << "Will scale over " << nthreads
#if defined(_OPENMP)
                << " OpenMP"
#else
                << " C++11"
#endif
                << " threads" << std::endl;
    }

    // count the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i) nelectron += atoms[i].atomic_number;
    //const auto ndocc = nelectron / 2;
    unsigned long naocc;
    if(nelectron % 2){
    	naocc = nelectron/2;
    }else{
    	double doublenaocc = static_cast <double>(nelectron)/2;
    	naocc = static_cast <unsigned long> (std::ceil(doublenaocc));
    }
    //unsigned long nbocc;
    auto nbocc = nelectron - naocc;

    cout << "# of electrons = " << nelectron << endl;

    // compute the nuclear repulsion energy
    auto enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij * xij + yij * yij + zij * zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    cout << "Nuclear repulsion energy = " << std::setprecision(15) << enuc
         << endl;

    libint2::Shell::do_enforce_unit_normalization(false);

    cout << "Atomic Cartesian coordinates (a.u.):" << endl;
    for (const auto& a : atoms)
      cout << a.atomic_number << " " << a.x << " " << a.y << " " << a.z
                << std::endl;

    BasisSet obs(basisname, atoms);
    cout << "orbital basis set rank = " << obs.nbf() << endl;

    // BASIS SET PROPERTIES QUANTUM NUMBERS NECESSARY FOR ANT.
    //vector<int> inputAOS = obs.shell2bf();
    vector<unsigned long> inputAOS = obs.shell2bf();
    //for (int i = inputAOS.size() - 1; i >= 0; i--)
    //    cout <<"inputAOS["<<i<<"]="<< inputAOS[i]<<"\n";
    for (int i = 0; i <= inputAOS.size() - 1; i++)
      cout <<"inputAOS["<<i<<"]="<< inputAOS[i]<<"\n";
    //vector<int> inputJAN = obs.shell2atom(atoms);
    vector<long> inputJAN = obs.shell2atom(atoms);
    //for (int i = inputJAN.size() - 1; i >= 0; i--)
    //        cout <<"inputJAN["<<i<<"]="<< inputJAN[i]<<"\n";
    for (int i = 0; i <= inputJAN.size() - 1; i++)
      cout <<"inputJAN["<<i<<"]="<< inputJAN[i]<<"\n";

    vector<size_t> inputShellT = obs.shell2bf();
    vector<size_t> inputShellC = obs.shell2bf();
    for (int i = 0; i < inputAOS.size()-1; i++){
      //inputShellT[i]=obs[i].am_symbol_to_l(obs[i].am_symbol());
      //inputShellT[i]=obs.max_l(obs[i].Shell());
      //inputShellT[i]=obs. (obs[i].Shell());
      inputShellT[i]=(inputAOS[i+1]-inputAOS[i]-1)/2;
      inputShellC[i]=0;
      cout <<"inputShellT["<<i<<"]="<< inputShellT[i]<<"\n";
    }

    inputShellT[inputAOS.size()-1]=(obs.nbf()-inputAOS[inputAOS.size()-1]-1)/2;
    inputShellC[inputAOS.size()-1]=0;
	cout <<"inputShellT["<<inputAOS.size()-1<<"]="<< inputShellT[inputAOS.size()-1]<<"\n";
	cout <<"obs.nbf()="<<obs.nbf()<<";inputAOS["<<inputAOS.size()-1<<"]="<< inputAOS[inputAOS.size()-1]<<"\n";

    //cout<<"Press ENTER to continue...";cin.get();
    // END BASIS SET PROPERTIES QUANTUM NUMBERS NECESSARY FOR ANT.

#ifdef HAVE_DENSITY_FITTING
    BasisSet dfbs;
    if (do_density_fitting) {
      dfbs = BasisSet(dfbasisname, atoms);
      cout << "density-fitting basis set rank = " << dfbs.nbf() << endl;
    }
#endif  // HAVE_DENSITY_FITTING

    /*** =========================== ***/
    /*** compute 1-e integrals       ***/
    /*** =========================== ***/

    // initializes the Libint integrals library ... now ready to compute
    libint2::initialize();

    // compute OBS non-negligible shell-pair list
    {
      obs_shellpair_list = compute_shellpair_list(obs);
      size_t nsp = 0;
      for (auto& sp : obs_shellpair_list) {
        nsp += sp.second.size();
      }
      cout << "# of {all,non-negligible} shell-pairs = {"
                << obs.size() * (obs.size() + 1) / 2 << "," << nsp << "}"
                << std::endl;
    }

    // compute one-body integrals
    auto S = compute_1body_ints<Operator::overlap>(obs)[0];
    auto T = compute_1body_ints<Operator::kinetic>(obs)[0];
    auto V = compute_1body_ints<Operator::nuclear>(obs, atoms)[0];
    Matrix H = T + V;
    //T.resize(0, 0); // COMMENTED ON 2018-01-23 TO BE PASSED TO initialize_Erkale_SCF.
    //V.resize(0, 0); // COMMENTED ON 2018-01-23 TO BE PASSED TO initialize_Erkale_SCF.

    strcpy(filename,basename);
    strcat(filename,".ham");
    writeEigenMatrixToCSVfile(filename, H);
    strcpy(filename,basename);
    strcat(filename,".over");
    writeEigenMatrixToCSVfile(filename, S);

    // compute orthogonalizer X such that X.transpose().S.X = I
    Matrix X, Xinv;
    double XtX_condition_number;  // condition number of "re-conditioned"
                                  // overlap obtained as Xinv.transpose().Xinv
    // one should think of columns of Xinv as the conditioned basis
    // Re: name ... cond # (Xinv.transpose().Xinv) = cond # (X.transpose() .
    // X)
    // by default assume can manage to compute with condition number of S <=
    // 1/eps
    // this is probably too optimistic, but in well-behaved cases even 10^11 is
    // OK
    double S_condition_number_threshold =
        1.0 / std::numeric_limits<double>::epsilon();
    std::tie(X, Xinv, XtX_condition_number) =
        conditioning_orthogonalizer(S, S_condition_number_threshold);

    cout << "XtX_condition_number = " << XtX_condition_number;
    cout << "S_condition_number_threshold = " << S_condition_number_threshold;

    Matrix Dalpha;
    Matrix Dbeta;
    Matrix D;
    Matrix Calpha_occ;
    Matrix Cbeta_occ;
    Matrix C_occ;
    Matrix evalsalpha;
    Matrix evalsbeta;
    Matrix evals;

    Matrix Calpha;
    Matrix Cbeta;
    Matrix XCalpha;
    Matrix XCbeta;
    Matrix Halpha;
    Matrix Hbeta;
    Matrix J;
    Matrix Kalpha;
    Matrix Kbeta;

    {  // use SOAD as the guess density
      const auto tstart = std::chrono::high_resolution_clock::now();

      auto D_minbs = compute_soad(atoms);  // compute guess in minimal basis
      //BasisSet minbs("STO-3G", atoms);
      BasisSet minbs(basisname, atoms);
      if (minbs == obs){
        Dalpha = D_minbs/2.0;
        Dbeta = D_minbs/2.0;
        D = D_minbs;
      }else {  // if basis != minimal basis, map non-representable SOAD guess
              // into the AO basis
              // by diagonalizing a Fock matrix
        cout << "projecting SOAD into AO basis ... ";
        auto Falpha = H;
        auto Fbeta = H;
        Falpha += compute_2body_fock_general(
            obs, D_minbs, minbs, true /* SOAD_D_is_shelldiagonal */,
            std::numeric_limits<double>::epsilon()  // this is cheap, no reason
                                                    // to be cheaper
            );
        Fbeta = Falpha;
        VectorOfMatrices F = {Falpha,Fbeta};

        // solve F C = e S C by (conditioned) transformation to F' C' = e C',
        // where
        // F' = X.transpose().F.X; the original C is obtained as C = X.C'
        Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_alpha(X.transpose() * Falpha * X);
        Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_beta(X.transpose() * Fbeta * X);
        auto Calpha = X * eig_solver_alpha.eigenvectors();
        auto Cbeta = X * eig_solver_beta.eigenvectors();

        // compute density, D = C(occ).C(occ)T
        Calpha_occ = Calpha.leftCols(naocc);
        Dalpha = Calpha_occ * Calpha_occ.transpose();
        Cbeta_occ = Cbeta.leftCols(nbocc);
        Dbeta = Cbeta_occ * Cbeta_occ.transpose();
        D = Dalpha + Dbeta;


        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;
        cout << "done (" << time_elapsed.count() << " s)" << std::endl;
      }
    }

    // pre-compute data for Schwarz bounds
    auto K = compute_schwarz_ints<>(obs);

// prepare for density fitting
#ifdef HAVE_DENSITY_FITTING
    std::unique_ptr<DFFockEngine> dffockengine(
        do_density_fitting ? new DFFockEngine(obs, dfbs) : nullptr);
#endif  // HAVE_DENSITY_FITTING

    cout << "*** =========================== ***\n";
    cout << "***          SCF loop           ***\n";
	cout << "*** =========================== ***\n";

    const auto maxiter = 10000;
//    const auto conv = 1e-12;
    const auto conv = 1e-8;
    auto iter = 0;
    auto rms_error = 1.0;
    double Falpha_err = 1.0;
    auto ediff_rel = 0.0;
    auto ehf = 0.0;
    //auto n2 = Dalpha.cols() * Dalpha.rows() + Dbeta.cols() * Dbeta.rows();
    auto n2 = Dalpha.cols() * Dalpha.rows();
    libint2::DIIS<Matrix> diis(2);  // start DIIS on second iteration

    // prepare for incremental Fock build ...
    Matrix D_diff = D;
    Matrix Dalpha_diff = Dalpha;
    Matrix Dbeta_diff = Dbeta;
    Matrix Falpha = H;
    Matrix Fbeta = H;
    VectorOfMatrices F = {Falpha,Fbeta};
    bool reset_incremental_fock_formation = false;
    bool incremental_Fbuild_started = false;
    double start_incremental_F_threshold = 1e-5;
    double next_reset_threshold = 0.0;
    size_t last_reset_iteration = 0;
    // ... unless doing DF, then use MO coefficients, hence not "incremental"
    if (do_density_fitting) start_incremental_F_threshold = 0.0;

    /*** =========================== ***/
    /***          DIIS loop          ***/
    /*** =========================== ***/
    bool do_diis = false;
    if(do_diis){
      do {
        const auto tstart = std::chrono::high_resolution_clock::now();
        ++iter;

        // Last iteration's energy and density
        auto ehf_last = ehf;
        Matrix D_last = D;
        Matrix Dalpha_last = Dalpha;
        Matrix Dbeta_last = Dbeta;

        if (not incremental_Fbuild_started &&
            rms_error < start_incremental_F_threshold) {
          incremental_Fbuild_started = true;
          reset_incremental_fock_formation = false;
          last_reset_iteration = iter - 1;
          next_reset_threshold = rms_error / 1e1;
          cout << "== started incremental fock build" << std::endl;
        }
        if (reset_incremental_fock_formation || not incremental_Fbuild_started) {
          Falpha = H;
          Fbeta = H;
          F = {Falpha,Fbeta};
          D_diff = D;
          Dalpha_diff = Dalpha;
          Dbeta_diff = Dbeta;
        }
        if (reset_incremental_fock_formation && incremental_Fbuild_started) {
          reset_incremental_fock_formation = false;
          last_reset_iteration = iter;
          next_reset_threshold = rms_error / 1e1;
          cout << "== reset incremental fock build" << std::endl;
        }

        // build a new Fock matrix
        if (not do_density_fitting) {
          // totally empirical precision variation, involves the condition number
          const auto precision_F = std::min(
              std::min(1e-3 / XtX_condition_number, 1e-7),
              std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon()));
          //F += compute_2body_fock_uhf(obs, Dalpha_diff, Dbeta_diff, precision_F, K);
          //F = std::plus<VectorOfMatrices>{}(F,compute_2body_fock_uhf(obs, Dalpha_diff, Dbeta_diff, precision_F, K));
          VectorOfMatrices pivotF = compute_2body_fock_uhf(obs, Dalpha_diff, Dbeta_diff, precision_F, K);
          Falpha = F[0] + pivotF[0];
          Fbeta = F[1] + pivotF[1];
        }
#if HAVE_DENSITY_FITTING
        else {  // do DF
          Falpha = H + dffockengine->compute_2body_fock_dfC(C_occ);
        }
#else
        else {
          assert(false);
        }  // do_density_fitting is true but HAVE_DENSITY_FITTING is not defined!
           // should not happen
#endif  // HAVE_DENSITY_FITTING

        //n2 = Dalpha.cols() * Dalpha.rows() + Dbeta.cols() * Dbeta.rows();
        // compute HF energy with the non-extrapolated Fock matrix
        //ehf = D.cwiseProduct(H + 0.5*(Falpha + Fbeta)).sum();
        ehf = Dalpha.cwiseProduct(H + Falpha).sum() + Dbeta.cwiseProduct(H + Fbeta).sum();
        ediff_rel = std::abs((ehf - ehf_last) / ehf);

        // compute SCF error
        D = Dalpha + Dbeta;
        Matrix FDalpha_comm = Falpha * Dalpha * S - S * Dalpha * Falpha ;
        Matrix FDbeta_comm = Fbeta * Dbeta * S - S * Dbeta * Fbeta;
        rms_error = 0.5*(FDalpha_comm.norm() + FDbeta_comm.norm()) / n2;
        if (rms_error < next_reset_threshold || iter - last_reset_iteration >= 8)
          reset_incremental_fock_formation = true;

        // DIIS extrapolate F
        Matrix Falpha_diis = Falpha;  // extrapolated F cannot be used in incremental Fock
                            // build; only used to produce the density
                            // make a copy of the unextrapolated matrix
        Matrix Fbeta_diis = Fbeta;  // extrapolated F cannot be used in incremental Fock
                            // build; only used to produce the density
                            // make a copy of the unextrapolated matrix
        diis.extrapolate(Falpha_diis, FDalpha_comm);
        diis.extrapolate(Fbeta_diis, FDbeta_comm);

        // solve F C = e S C by (conditioned) transformation to F' C' = e C',
        // where
        // F' = X.transpose().F.X; the original C is obtained as C = X.C'
        Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_alpha(X.transpose() * Falpha_diis *
                                                       X);
        Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_beta(X.transpose() * Fbeta_diis *
                                                       X);
        evalsalpha = eig_solver_alpha.eigenvalues();
        auto Calpha = X * eig_solver_alpha.eigenvectors();
        evalsbeta = eig_solver_beta.eigenvalues();
        auto Cbeta = X * eig_solver_beta.eigenvectors();

        // compute density, D = C(occ).C(occ)T
        Calpha_occ = Calpha.leftCols(naocc);
        Dalpha = Calpha_occ * Calpha_occ.transpose();
        Dalpha_diff = Dalpha - Dalpha_last;
        Cbeta_occ = Cbeta.leftCols(nbocc);
        Dbeta = Cbeta_occ * Cbeta_occ.transpose();
        Dbeta_diff = Dbeta - Dbeta_last;

        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        if (iter == 1)
          cout << "\nIter         E(HF)                 D(E)/E         "
                       "RMS([F,D])/nn       Time(s)\n";
        printf(" %02d %20.12f %20.12e %20.12e %10.5lf\n", iter, ehf + enuc,
               ediff_rel, rms_error, time_elapsed.count());

      } while (((ediff_rel > conv) || (rms_error > conv)) && (iter < maxiter));
    }
    /*** =========================== ***/
    /***        END DIIS loop        ***/
    /*** =========================== ***/

    cout << "*** =========================== ***\n";
    cout << "***          ANT loop           ***\n";
	cout << "*** =========================== ***\n";
	cout << "Initialize some parameters to be passed to ANTC\n";

	bool UHF = true;
	//char *inputjobname;
	//strcpy(inputjobname,filename);
	cout << "basename = "<<basename<<"\n";
	char inputjobname[strlen(basename)];
	char erkalejobname[strlen(basename)];
	//inputjobname << filename;
	//for( size_t i = 0; i < strlen(basename); i++ ) // assign each element of nameInput to TESCStudent.name explicitly
	//{
	//  inputjobname[ i ] = basename[ i ];
	//}

	// ORIGINAL CODE, WORKED BEFORE 2018-12-31.
	//strcpy(inputjobname,basename); // ORIGINAL CODE, WORKED BEFORE 2018-12-31.
	strcpy(inputjobname,fullname);
	strcpy(erkalejobname,basename);
	strcat(erkalejobname,"_");
	strcat(erkalejobname,basisname);

	cout << "inputjobname = "<<inputjobname<<"\n";
	cout << "erkalejobname = "<<erkalejobname<<"\n";
	unsigned long inputjobname_len = strlen(inputjobname);
	unsigned long inputJCycle; // Later defined as iter-1 for every cycle.
	unsigned long inputNCycles = maxiter;

    unsigned long NBasis = obs.nbf();
    unsigned long inputNSpin;
    unsigned long inputNAtoms = atoms.size();
    unsigned long inputNShell = inputAOS.size();

    unsigned long inputNAE;
    unsigned long inputNBE;
    unsigned long inputNE;
    int mult;

    if(nelectron % 2){
    	inputNAE = nelectron/2;
    	inputNSpin = 1;
    }else{
    	double doubleNAE = static_cast <double>(nelectron)/2;
    	inputNAE = static_cast <unsigned long> (std::ceil(doubleNAE));
    	inputNSpin = 2;
    }
    inputNBE = nelectron - inputNAE;
    inputNE = nelectron;
    mult = inputNSpin;


    cout << "BEFORE initialize_Erkale_SCF(basename,NBasis,nelectron,mult,S,T,V,H);\n";
    //ErkaleBasisSet * erkalebasis;
    ErkaleBasisSet erkalebasis;
    //ErkaleBasisSet *erkalebasis = new ErkaleBasisSet;
    Settings erkalesettings;
    //initialize_Erkale_SCF(string basename, unsigned long int Nbf, int Nel, int mult, Matrix& S, Matrix& T, Matrix& Vnuc, Matrix& Hcore )
    SCF ErkaleSCF = initialize_Erkale_SCF(erkalejobname,erkalebasis,erkalesettings,NBasis,nelectron,mult,S,T,V,H);
    cout << "\n******************************************************\n";
    cout << "\nPRINT ERKALEBASIS in runant\n";
    cout << "\n******************************************************\n";
    erkalebasis.print(true);
    cout << "\n******************************************************\n";
    cout << "\nEND PRINT ERKALEBASIS in runant\n";
    cout << "\n******************************************************\n";
    /*
    cout << "\n******************************************************\n";
    cout << "\nPRINT *ERKALEBASIS in runant\n";
    cout << "\n******************************************************\n";
    *erkalebasis->print(true);
    cout << "\n******************************************************\n";
    cout << "\nEND PRINT *ERKALEBASIS in runant\n";
    cout << "\n******************************************************\n";
    */
    cout << "AFTER initialize_Erkale_SCF(basename,NBasis,nelectron,mult,S,T,V,H);\n";




    Matrix FockHWA = Matrix::Zero(NBasis,NBasis);
    Matrix FockHWB = Matrix::Zero(NBasis,NBasis);
    Matrix GibbsYA = Matrix::Zero(NBasis,NBasis);
    Matrix GibbsYB = Matrix::Zero(NBasis,NBasis);
    Matrix GibbsYKernel1A = Matrix::Zero(NBasis,NBasis);
    Matrix GibbsYKernel1B = Matrix::Zero(NBasis,NBasis);
    Matrix GibbsYKernel2A = Matrix::Zero(NBasis,NBasis);
    Matrix GibbsYKernel2B = Matrix::Zero(NBasis,NBasis);

    cout << "BEFORE DECLARE inputIAN, inputAtmChg, inputAtmCo.\n";
    vector<unsigned long> inputIAN;
    inputIAN.resize(atoms.size());
    //double inputAtmCo[3][atoms.size()]={0.0};
    vector<double> inputAtmChg;
    inputAtmChg.resize(atoms.size());

    //double** inputAtmCo = new double*[3]; // rowCount = 3
    //for(int i = 0; i < 3; ++i)
    //  inputAtmCo[i] = new double[atoms.size()]; // colCount = atoms.size()

    double **inputAtmCo;
    inputAtmCo = alloc_2d_double(atoms.size(),3);
    //double *inputAtmCo = new double[3*inputAtmCo];
    //for(int i = 0; i < 3; ++i)
    //  inputAtmCo[i] = new double[inputAtmCo]; // colCount = atoms.size()

    //vector<vector<double> > inputAtmCo;
    //inputAtmCo.resize(3);
    //for (int i = 0; i < 3; ++i)
    //   inputAtmCo[i].resize(atoms.size());

    //Matrix inputAtmCo = Matrix::Zero(3,atoms.size());
    //inputAtmCo[0].data();
    //inputAtmCo.data()[0];
    cout << "BEFORE FILL inputIAN, inputAtmChg, inputAtmCo.\n";
    //double inputAtmChg[atoms.size()];
    //int iat = 0;
    for(int iat = 0; iat < atoms.size();++iat){
      cout << "atoms["<<iat<<"].atomic_number = "<<atoms[iat].atomic_number<<"\n";
      inputIAN[iat]=atoms[iat].atomic_number;
      inputAtmCo[iat][0]=atoms[iat].x;
      inputAtmCo[iat][1]=atoms[iat].y;
      inputAtmCo[iat][2]=atoms[iat].z;
      inputAtmChg[iat]=atoms[iat].atomic_number;
      cout << "iat,inputIAN,inputAtmChg : "<<iat<<","<<inputIAN[iat]<<","<<inputAtmChg[iat]<<"\n";
      cout << "inputAtmCo : ("<<inputAtmCo[iat][0]<<","<<inputAtmCo[iat][1]<<","<<inputAtmCo[iat][2]<<")\n";
    }
    cout << "AFTER DECLARE inputIAN, inputAtmChg, inputAtmCo.\n";
    cout << "BEFORE DECLARE Darray.\n";
    //auto Darray = D;
    //auto DA = Dalpha.data();
    //auto DB = Dbeta.data();
    double * DA = Dalpha.data();
    double * DB = Dbeta.data();
    auto pivHam = H.data();
    auto pivFockA = Falpha.data();
    auto pivFockB = Fbeta.data();
    auto pivCoulomb = Falpha.data();
    auto pivExchangeA = Falpha.data();
    auto pivExchangeB = Fbeta.data();
    auto pivOver = S.data();

    auto outFockHWA = FockHWA.data();
    auto outFockHWB = FockHWB.data();
    auto outGibbsYA = GibbsYA.data();
    auto outGibbsYB = GibbsYB.data();
    auto outGibbsYKernel1A = GibbsYKernel1A.data();
    auto outGibbsYKernel1B = GibbsYKernel1B.data();
    auto outGibbsYKernel2A = GibbsYKernel2A.data();
    auto outGibbsYKernel2B = GibbsYKernel2B.data();
    //auto pivOver = compute_1body_ints<Operator::overlap>(obs)[0];
    cout << "AFTER DECLARE pivOver.\n";
    auto denerrj = 1e-5;
    auto Crit = 1e-5;
    bool ANTOn = true;
    cout << "End initialize some parameters to be passed to ANTC\n";

    bool do_ant = true;
    if(do_ant){
    	cout << "ENTER do_ant IF\n";
        do {
        	cout << "ENTER do_ant DO\n";
          const auto tstart = std::chrono::high_resolution_clock::now();
          ++iter;
          inputJCycle = iter-1;

          // Last iteration's energy and density
          auto ehf_last = ehf;
          Matrix Dalpha_last = Dalpha;
          Matrix Dbeta_last = Dbeta;
          Matrix D_last = D;
          auto Falpha_last = Falpha;
          auto Fbeta_last = Fbeta;

          if (not incremental_Fbuild_started &&
              rms_error < start_incremental_F_threshold) {
            incremental_Fbuild_started = true;
            reset_incremental_fock_formation = false;
            last_reset_iteration = iter - 1;
            next_reset_threshold = rms_error / 1e1;
            if (debugflag)cout << "== started incremental fock build" << std::endl;
            if (debugflag)cout << "(not incremental_Fbuild_started && rms_error < start_incremental_F_threshold)\n";
            //cout<<"Press ENTER to continue...";cin.get();
          }
          if (debugflag)cout << "AFTER if (not incremental_Fbuild_started &&"<<
        		  	   "rms_error < start_incremental_F_threshold)\n";
          if (reset_incremental_fock_formation || not incremental_Fbuild_started) {
            Falpha = H;
            Fbeta = H;
            F = {Falpha,Fbeta};
            Dalpha_diff = Dalpha;
            Dbeta_diff = Dbeta;
            D_diff = D;
            if (debugflag)cout<<"(reset_incremental_fock_formation || not incremental_Fbuild_started)\n";
            //cout<<"Press ENTER to continue...";cin.get();
          }
          if (debugflag)cout << "AFTER if (reset_incremental_fock_formation || not incremental_Fbuild_started)\n";
          if (reset_incremental_fock_formation && incremental_Fbuild_started) {
            reset_incremental_fock_formation = false;
            last_reset_iteration = iter;
            next_reset_threshold = rms_error / 1e1;
            if (debugflag)cout << "== reset incremental fock build" << std::endl;
            //cout<<"Press ENTER to continue...";cin.get();
          }
          if (debugflag)cout << "AFTER if (reset_incremental_fock_formation && incremental_Fbuild_started)\n";

          // build a new Fock matrix
          if (not do_density_fitting) {
        	  if (debugflag)cout << "INSIDE if (not do_density_fitting)\n";
            // totally empirical precision variation, involves the condition number
            const auto precision_F = std::min(
                std::min(1e-3 / XtX_condition_number, 1e-7),
                std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon()));

            if (debugflag){
            cout << "rms_error = " << rms_error<<"\n";
            cout << "1e-3 / XtX_condition_number = " << 1e-3 / XtX_condition_number<<"\n";
            cout << "XtX_condition_number = " << XtX_condition_number<<"\n";
            cout << "std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon()) = " << std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon())<<"\n";
            cout << "INSIDE if (not do_density_fitting) BEFORE compute_2body_fock\n";
            cout << "START printEigenmatrix(D_diff)\n";
            //printEigenmatrix(D_diff);
            cout << "END printEigenmatrix(D_diff)\n";
            cout << "START printEigenmatrix(K)\n";
            //printEigenmatrix(K);
            cout << "END printEigenmatrix(K)\n";
            cout << "START printEigenmatrix(F)\n";
            //printEigenmatrix(F);
            cout << "END printEigenmatrix(F)\n";
            cout << "precision_F = " << precision_F << "\n";
            }

            //F += compute_2body_fock(obs, D_diff, precision_F, K);
            auto F_add = compute_2body_fock_uhf(obs, Dalpha_diff, Dbeta_diff, precision_F, K);
            if (debugflag){
            cout << "START printEigenmatrix(F_add)\n";
            //printEigenmatrix(F_add);
            cout << "END printEigenmatrix(F_add)\n";
            }
            //F = F + F_add;
            //F = std::plus<VectorOfMatrices*>{}(F,F_add);
            //Falpha = 0.9*Falpha_last + 0.1*(F[0] + F_add[0]);
            //Fbeta = 0.9*Fbeta_last + 0.1*(F[1] + F_add[1]);
            Falpha = F[0] + F_add[0];
            Fbeta = F[1] + F_add[1];
            F = {Falpha,Fbeta};
            if (debugflag){
            cout << "START printEigenmatrix(F)\n";
            //printEigenmatrix(F);
            cout << "END printEigenmatrix(F)\n";
            }




          }
    #if HAVE_DENSITY_FITTING
          else {  // do DF
        	  if (debugflag)cout << "INSIDE HAVE_DENSITY_FITTING\n";
              //F = H + dffockengine->compute_2body_fock_dfC(C_occ);
        	  Falpha = H + dffockengine->compute_2body_fock_dfC(Calpha_occ);
        	  Fbeta = H + dffockengine->compute_2body_fock_dfC(Cbeta_occ);
        	  F = {Falpha,Fbeta};
        	  //cout<<"Press ENTER to continue...";cin.get();
          }
    #else
          else {
        	  if (debugflag)cout << "INSIDE NOT HAVE_DENSITY_FITTING\nCOMMENTED assert(false);\n";
            //assert(false);

        	  // VARIABLES QUE ES NECESARIO ACTUALIZAR.
        	  cout<<"VARIABLES QUE ES NECESARIO ACTUALIZAR.\n";

        	  const auto precision_F = std::min(
        			  std::min(1e-3 / XtX_condition_number, 1e-7),
					  std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon()));

        	  if (debugflag){
        	  	cout << "rms_error = " << rms_error<<"\n";
        	    cout << "1e-3 / XtX_condition_number = " << 1e-3 / XtX_condition_number<<"\n";
        	    cout << "XtX_condition_number = " << XtX_condition_number<<"\n";
        	    cout << "std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon()) = " << std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon())<<"\n";
        	    cout << "INSIDE if (not do_density_fitting) BEFORE compute_2body_fock\n";
        	    cout << "START printEigenmatrix(D_diff)\n";
        	    //printEigenmatrix(D_diff);
        	    cout << "END printEigenmatrix(D_diff)\n";
        	    cout << "START printEigenmatrix(K)\n";
        	    //printEigenmatrix(K);
        	    cout << "END printEigenmatrix(K)\n";
        	    cout << "START printEigenmatrix(F)\n";
        	    //printEigenmatrix(F);
        	    cout << "END printEigenmatrix(F)\n";
        	    cout << "precision_F = " << precision_F << "\n";
        	  }

        	  Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_alpha(X.transpose() * Falpha * X);
        	  Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_beta(X.transpose() * Fbeta * X);
        	  //auto Calpha = X * eig_solver_alpha.eigenvectors();
        	  //auto Cbeta = X * eig_solver_beta.eigenvectors();
        	  Matrix Calpha = X * eig_solver_alpha.eigenvectors();
        	  Matrix Cbeta = X * eig_solver_beta.eigenvectors();
        	  evalsalpha = eig_solver_alpha.eigenvalues();
        	  evalsbeta = eig_solver_beta.eigenvalues();

        	  // compute density, D = C(occ).C(occ)T
        	  Calpha_occ = Calpha.leftCols(naocc);
        	  //Dalpha = Calpha_occ * Calpha_occ.transpose();
        	  Cbeta_occ = Cbeta.leftCols(nbocc);
        	  //Dbeta = Cbeta_occ * Cbeta_occ.transpose();
        	  //D = Dalpha + Dbeta;

        	  //F += compute_2body_fock(obs, D_diff, precision_F, K);
        	  auto JK_add = compute_2body_JK_uhf(obs, Dalpha_diff, Dbeta_diff, precision_F, K);
        	  if (debugflag){
        	  	cout << "START printEigenmatrix(JK_add)\n";
        	    //printEigenmatrix(F_add);
        	    cout << "END printEigenmatrix(JK_add)\n";
        	  }

        	  auto J = JK_add[0];
        	  auto Kalpha = JK_add[1];
        	  auto Kbeta = JK_add[2];

        	  Matrix XCalpha = Matrix::Zero(NBasis,NBasis);
        	  Matrix XCbeta = Matrix::Zero(NBasis,NBasis);

        	  if (debugflag)cout << "CALLING DFT FROM ERKALE\n";

        	  cout << "START printEigenmatrix(evalsalpha)\n";
        	  printEigenmatrix(evalsalpha);
        	  cout << "END printEigenmatrix(evalsalpha)\n";

              //double Ea = 0.0;
              //double Eb = 0.0;
        	  //std::vector<double> Ea(evalsalpha.diagonal().data(),evalsalpha.diagonal().data() + evalsalpha.diagonal().size());
        	  //std::vector<double> Eb(evalsbeta.diagonal().data(),evalsbeta.diagonal().data() + evalsbeta.diagonal().size());
        	  // evalsalpha y evalsbeta NO SON MATRICES DIAGONALES, SON VECTORES.
        	  std::vector<double> Ea(evalsalpha.data(),evalsalpha.data() + evalsalpha.size());
        	  std::vector<double> Eb(evalsbeta.data(),evalsbeta.data() + evalsbeta.size());

        	  cout <<"evalsalpha.size() = "<<evalsalpha.size()<<"\n";

              cout <<"Ea";
              printStdVector(Ea);
              cout <<"Eb";
              printStdVector(Eb);



              cout << "START printEigenmatrixCleanFmt(Falpha)\n";
              printEigenmatrixCleanFmt(Falpha);
              cout << "END printEigenmatrixCleanFmt(Falpha)\n";

              arma::mat Fa = example_cast_arma(Falpha);
              Fa.print(std::cout);
              cout << "\n";

              cout << "START printEigenmatrixCleanFmt(Fbeta)\n";
              printEigenmatrixCleanFmt(Fbeta);
              cout << "END printEigenmatrixCleanFmt(Fbeta)\n";

              arma::mat Fb = example_cast_arma(Fbeta);
              Fb.print(std::cout);
              cout << "\n";

              cout<<"END VARIABLES QUE ES NECESARIO ACTUALIZAR.\n";

              cout << "BEFORE wrap_Fock_UDFT(ErkaleSCF...);\n";
              //unsigned long int Nbf = obs.nbf(); // or NBasis
              //wrap_Fock_UDFT(SCF ErkaleSCF, string basename, const double Ea, const double Eb, const Matrix& Calpha, const Matrix& Cbeta, const Matrix& Halpha, const Matrix& Hbeta, const Matrix& Dalpha, const Matrix& Dbeta, const Matrix& J, const Matrix& Kalpha, const Matrix& Kbeta, const Matrix& XCalpha, const Matrix& XCbeta)
              wrap_Fock_UDFT(ErkaleSCF,erkalebasis,erkalesettings,erkalejobname,NBasis,Ea,Eb,naocc,nbocc,S,Calpha,Cbeta,Calpha_occ,Cbeta_occ,Falpha,Fbeta,Dalpha,Dbeta,J,Kalpha,Kbeta,XCalpha,XCbeta);
              cout << "AFTER wrap_Fock_UDFT(ErkaleSCF...);\n";
              //erkalebasis.get_Nnuc();
              //ErkaleSCF.


              Falpha = F[0] + J + Kalpha + XCalpha;
              Fbeta = F[1] + J + Kbeta + XCbeta;
              F = {Falpha,Fbeta};
              if (debugflag){
              	cout << "START printEigenmatrix(F)\n";
                //printEigenmatrix(F);
                cout << "END printEigenmatrix(F)\n";
               }

          }  // do_density_fitting is true but HAVE_DENSITY_FITTING is not defined!
             // should not happen
    #endif  // HAVE_DENSITY_FITTING

          cout << "AFTER if (not do_density_fitting)\n";

          auto Falpha_diff = (Falpha - Falpha_last).cwiseAbs();
          Eigen::MatrixXf::Index maxRow, maxCol;
          Falpha_err = Falpha_diff.maxCoeff(&maxRow, &maxCol);

          // MOVE THIS AFTER ANTC CALL TO EVALUATE PROPER ehf AND ediff_rell.
          // compute HF energy with the non-extrapolated Fock matrix
          //ehf = Dalpha.cwiseProduct(H + Falpha).sum() + Dbeta.cwiseProduct(H + Fbeta).sum();
          ehf = Dalpha.cwiseProduct(H + Falpha).sum();
          ediff_rel = std::abs((ehf - ehf_last) / ehf);

          if (debugflag)cout << "BEFORE compute SCF error\n";
          // compute SCF error
          D = Dalpha + Dbeta;
          //Matrix FD_comm = F * D * S - S * D * F;
          Matrix FDalpha_comm = Falpha * Dalpha * S - S * Dalpha * Falpha;
          Matrix FDbeta_comm = Fbeta * Dbeta * S - S * Dbeta * Fbeta;
          //rms_error = (FDalpha_comm.norm() + FDbeta_comm.norm())/ n2;
          rms_error = FDalpha_comm.norm() / n2;
          //rms_error = 0.5*(FDalpha_comm.norm() + FDbeta_comm.norm()) / n2;

          if (rms_error < next_reset_threshold || iter - last_reset_iteration >= 8)
            reset_incremental_fock_formation = true;

          cout << "BEFORE DIIS extrapolate F\n";
          // DIIS extrapolate F
          Matrix Falpha_diis = Falpha;  // extrapolated F cannot be used in incremental Fock
                              // build; only used to produce the density
                              // make a copy of the unextrapolated matrix
          Matrix Fbeta_diis = Fbeta;
          diis.extrapolate(Falpha_diis, FDalpha_comm);
          diis.extrapolate(Fbeta_diis, FDbeta_comm);

          cout << "BEFORE do_diis REPLACED BY THE ANT EXECUTION.\n";
          // THIS PART IS REPLACED BY THE ANT EXECUTION.
          do_diis = false;
          if(do_diis){
        	  cout << "INSIDE do_diis SHOULD NOT BE HERE BECAUSE REPLACED BY THE ANT EXECUTION.\n";
            // solve F C = e S C by (conditioned) transformation to F' C' = e C',
            // where
            // F' = X.transpose().F.X; the original C is obtained as C = X.C'
            Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_alpha(X.transpose() * Falpha_diis *
                                                             X);
            Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_beta(X.transpose() * Fbeta_diis *
                                                             X);
            evalsalpha = eig_solver_alpha.eigenvalues();
            auto Calpha = X * eig_solver_alpha.eigenvectors();
            evalsbeta = eig_solver_beta.eigenvalues();
            auto Cbeta = X * eig_solver_beta.eigenvectors();

            // compute density, D = C(occ).C(occ)T
            Calpha_occ = Calpha.leftCols(naocc);
            Dalpha = Calpha_occ * Calpha_occ.transpose();
            Dalpha_diff = Dalpha - Dalpha_last;
            Cbeta_occ = Cbeta.leftCols(nbocc);
            Dbeta = Cbeta_occ * Cbeta_occ.transpose();
            Dbeta_diff = Dbeta - Dbeta_last;
            D_diff = Dalpha_diff + Dbeta_diff;
            cout<<"INSIDE DO_DIIS. Press ENTER to continue...";cin.get();
          }else{
        	  if (debugflag)cout << "INSIDE  NOT do_diis BECAUSE REPLACED BY THE ANT EXECUTION.\n";
        	  Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_alpha(X.transpose() * Falpha_diis *
        	                                                               X);
        	  Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_beta(X.transpose() * Fbeta_diis *
        	          	                                                   X);
        	  evalsalpha = eig_solver_alpha.eigenvalues();
        	  auto Calpha = X * eig_solver_alpha.eigenvectors();
        	  Calpha_occ = Calpha.leftCols(naocc);
        	  evalsbeta = eig_solver_beta.eigenvalues();
        	  auto Cbeta = X * eig_solver_beta.eigenvectors();
        	  Cbeta_occ = Cbeta.leftCols(nbocc);


        	  if (debugflag)cout << "BEFORE DECLARE Darray.\n";
        	  //auto Darray = D;
        	  DA = Dalpha.data();

        	  if (debugflag)cout << "PRINT ARRAY DA BEFORE ANTC.\n";
        	  printDoubleMatrix(DA,NBasis,NBasis);


        	  DB = Dbeta.data();

        	  if (debugflag)cout << "PRINT ARRAY DB BEFORE ANTC.\n";
        	  printDoubleMatrix(DB,NBasis,NBasis);


        	  pivHam = H.data();
        	  pivFockA = Falpha.data();
        	  pivFockB = Fbeta.data();
        	  pivCoulomb = Falpha.data();
        	  pivExchangeA = Falpha.data();
        	  pivExchangeB = Fbeta.data();
        	  pivOver = S.data();
        	  //auto pivOver = compute_1body_ints<Operator::overlap>(obs)[0];
        	  if (debugflag)cout << "AFTER DECLARE pivOver.\n";

        	  cout << "BEFORE CALL antc\n";
        	  //cout<<"Press ENTER to continue...";cin.get();

        	  //obs.
//        	  HartreeFock::matrix_multiply(&A[0], &rows_A, &cols_A, &B[0], &rows_B, &cols_B, &C[0], &rows_C, &cols_C);
//        	  HartreeFock::antc(bool *UHF, int *JCycle, int *inputNCycles,char *inputjobname, int *inputjobname_len,double *D,double *pivHam,double *pivFock,double *pivCoulomb,double *pivExchange,double *pivOver,
//        	       int *NBasis,int *inputNSpin,int *inputNAtoms,int *inputNShell,int *inputJAN,int *inputAOS,int *inputShellT,int *inputShellC,
//        	       double *inputNAE,double *inputNBE,double *inputNE,int *inputIAN,double *inputAtmCo,double *inputAtmChg,
//        	       int *IRwH,int *IRwPA,int *IRwPB,int *IRwFA,int *IRwFB,int *IRwS1,int *IRwEig,double *denerrj,double *Crit,bool *ANTOn);
        	  /*
        	  HartreeFock::antc(&UHF, &inputJCycle, &inputNCycles,&inputjobname[0],&inputjobname_len,
        			  &Darray[0],
        			  &pivHam[0],
					  &pivFock[0],
					  &pivCoulomb[0],
					  &pivExchange[0],
					  &pivOver[0],
					  &NBasis,&inputNSpin,&inputNAtoms,&inputNShell,
					  &inputJAN[0],
					  &inputAOS[0],
					  &inputShellT[0],
					  &inputShellC[0],
					  &inputNAE,&inputNBE,&inputNE,
					  &inputIAN[0],
					  &inputAtmCo[0][0],
					  &inputAtmChg[0],
					  &denerrj,&Crit,&ANTOn);
			  */
        	  HartreeFock::antc(&UHF,&inputJCycle,&inputNCycles,&inputjobname[0],&inputjobname_len,
              		  &NBasis,&inputNSpin,&inputNAtoms,&inputNShell,
					  &DA[0],
					  &DB[0],
					  &pivHam[0],
					  &pivFockA[0],
					  &pivFockB[0],
					  &pivCoulomb[0],
					  &pivExchangeA[0],
					  &pivExchangeB[0],
					  &pivOver[0],
					  &outFockHWA[0],
					  &outFockHWB[0],
					  &outGibbsYA[0],
					  &outGibbsYB[0],
					  &outGibbsYKernel1A[0],
					  &outGibbsYKernel1B[0],
					  &outGibbsYKernel2A[0],
					  &outGibbsYKernel2B[0],
					  &inputIAN[0],
					  &inputJAN[0],
					  &inputAOS[0],
					  &inputShellT[0],
					  &inputShellC[0],
					  &inputNAE,&inputNBE,&inputNE,
					  //&inputAtmCo.data()[0][0],
					  &(inputAtmCo[0][0]),
					  &inputAtmChg[0],
					  &denerrj,&Crit,&ANTOn);


        	  if (debugflag){
        	  cout << "START printDoubleMatrix(DA)\n";
        	  cout << "NBasis = "<<NBasis<<"\n";
        	  printDoubleMatrix(DA,static_cast<double>(NBasis),static_cast<double>(NBasis));
        	  //print_matrix(&DA[0],NBasis,NBasis);
        	  cout << "END printEigenmatrixOctaveFmt(DA)\n";
        	  cout << "START printDoubleMatrix(DB)\n";
        	  cout << "NBasis = "<<NBasis<<"\n";
        	  printDoubleMatrix(DB,NBasis,NBasis);
        	  //print_matrix(&DB[0],NBasis,NBasis);
        	  cout << "END printDoubleMatrix(DB)\n";
        	  }


        	  Dalpha = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(DA,NBasis,NBasis);
        	  Dalpha_diff = Dalpha - Dalpha_last;
        	  Dbeta = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(DB,NBasis,NBasis);
        	  Dbeta_diff = Dbeta - Dbeta_last;
        	  D = Dalpha + Dbeta;
        	  D_diff = Dalpha_diff + Dbeta_diff;


        	  cout << "START printEigenmatrixCleanFmt(Dalpha)\n";
        	  printEigenmatrixCleanFmt(Dalpha);
        	  cout << "END printEigenmatrixCleanFmt(Dalpha)\n";
        	  cout << "START printEigenmatrixCleanFmt(Dbeta)\n";
        	  printEigenmatrixCleanFmt(Dbeta);
        	  cout << "END printEigenmatrixCleanFmt(Dbeta)\n";
        	  // MOVED HERE TO EVALUATE PROPER ehf AND ediff_rel WITH UPDATED DENSITIES.
        	  // compute HF energy with the non-extrapolated Fock matrix
        	  /*
        	  ehf = Dalpha.cwiseProduct(H + Falpha).sum() + Dbeta.cwiseProduct(H + Fbeta).sum();
        	  ediff_rel = std::abs((ehf - ehf_last) / ehf);
        	  cout <<"(((ediff_rel > conv) || (rms_error > conv)) && (iter < maxiter))\n";
        	  cout <<"((("<<ediff_rel<<">"<<conv<<") || ("<<rms_error<<">"<<conv<<")) && ("<<iter<<" < "<<maxiter<<"))\n";
        	  */


        	  cout << "AFTER CALL antc\n";

        	  //cout<<"Press ENTER to continue...";cin.get();

          }

          // END OF PART TO BE REPLACED BY THE ANT EXECUTION.

          // MOVED HERE TO EVALUATE PROPER ehf AND ediff_rel WITH UPDATED DENSITIES.
          // compute HF energy with the non-extrapolated Fock matrix
          /*
          ehf = Dalpha.cwiseProduct(H + Falpha).sum() + Dbeta.cwiseProduct(H + Fbeta).sum();
          ediff_rel = std::abs((ehf - ehf_last) / ehf);
          */

          const auto tstop = std::chrono::high_resolution_clock::now();
          const std::chrono::duration<double> time_elapsed = tstop - tstart;

          if (iter == 1)
            cout << "\nIter         E(HF)                 D(E)/E         "
                         "RMS([F,D])/nn       Falpha_err       Time(s)\n";
          cout << "-------------------------------------------------------------------\n";
          cout << "------------------- DISPLAY SCF ANT CONVERGENCE -------------------\n";
          cout << "-------------------------------------------------------------------\n";
          cout << "\nIter         E(HF)                 D(E)/E         "
                       "RMS([F,D])/nn       Falpha_err       Time(s)\n";
          printf(" %02d %20.12f %20.12e %20.12e %20.12e %10.5lf\n", iter, ehf + enuc,
                 ediff_rel, rms_error, Falpha_err, time_elapsed.count());

          //cout <<"(((ediff_rel > conv) || (rms_error > conv)) && (iter < maxiter))\n";
          //cout <<"((("<<ediff_rel<<">"<<conv<<") || ("<<rms_error<<">"<<conv<<")) && ("<<iter<<" < "<<maxiter<<"))\n";
          cout <<"(((ediff_rel > conv) || (Falpha_err > conv)) && (iter < maxiter))\n";
          cout <<"((("<<ediff_rel<<">"<<conv<<") || ("<<Falpha_err<<">"<<conv<<")) && ("<<iter<<" < "<<maxiter<<"))\n";
          cout << "-------------------------------------------------------------------\n";
          cout << "----------------- END DISPLAY SCF ANT CONVERGENCE -----------------\n";
          cout << "-------------------------------------------------------------------\n";

          cout << "EXIT do_ant DO\n";

          //if((ediff_rel <= conv) && (rms_error <= conv)){
          if((ediff_rel <= conv) && (Falpha_err <= conv) && (iter >= 20)){
        	  //iter = maxiter - 1;
        	  //inputNCycles = iter; // WITH THIS LINE THE LAST CYCLE EVALUATION WILL BE FORCED IN NEXT CALL.
        	  inputJCycle++; // WITH THIS LINE THE LAST CYCLE EVALUATION WILL BE FORCED.
        	  inputNCycles = inputJCycle; // WITH THIS LINE THE LAST CYCLE EVALUATION WILL BE FORCED.
        	  cout << "---------------------------------------------------------\n";
        	  cout << "-------------------- LAST ANTC CYCLE --------------------\n";
        	  cout << "---------------------------------------------------------\n";




        	  HartreeFock::antc(&UHF,&inputJCycle,&inputNCycles,&inputjobname[0],&inputjobname_len,
        			  &NBasis,&inputNSpin,&inputNAtoms,&inputNShell,
        	  		  &DA[0],
        	  		  &DB[0],
        	  		  &pivHam[0],
					  &pivFockA[0],
        	  		  &pivFockB[0],
        	  		  &pivCoulomb[0],
        	  		  &pivExchangeA[0],
        	  		  &pivExchangeB[0],
        	  		  &pivOver[0],
					  &outFockHWA[0],
					  &outFockHWB[0],
					  &outGibbsYA[0],
					  &outGibbsYB[0],
					  &outGibbsYKernel1A[0],
					  &outGibbsYKernel1B[0],
					  &outGibbsYKernel2A[0],
					  &outGibbsYKernel2B[0],
        	  		  &inputIAN[0],
        	  		  &inputJAN[0],
        	  		  &inputAOS[0],
        	  		  &inputShellT[0],
        	  		  &inputShellC[0],
        	  		  &inputNAE,&inputNBE,&inputNE,
        	  		  //&inputAtmCo.data()[0][0],
        	  		  &(inputAtmCo[0][0]),
        	  		  &inputAtmChg[0],
        	  		  &denerrj,&Crit,&ANTOn);



        	  printf("---------------------------------------------------------\n");
        	  printf("----------------- AFTER LAST ANTC CALL ------------------\n");
        	  printf("---------------------------------------------------------\n");

        	  cout.flush();

        	  Dalpha = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(DA,NBasis,NBasis);
        	  Dalpha_diff = Dalpha - Dalpha_last;
        	  Dbeta = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(DB,NBasis,NBasis);
        	  Dbeta_diff = Dbeta - Dbeta_last;
        	  D = Dalpha + Dbeta;
        	  D_diff = Dalpha_diff + Dbeta_diff;

        	  cout << "START printEigenmatrixCleanFmt(Dalpha)\n";
        	  printEigenmatrixCleanFmt(Dalpha);
        	  cout << "END printEigenmatrixCleanFmt(Dalpha)\n";
        	  cout << "START printEigenmatrixCleanFmt(Dbeta)\n";
        	  printEigenmatrixCleanFmt(Dbeta);
        	  cout << "END printEigenmatrixCleanFmt(Dbeta)\n";

        	  FockHWA = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outFockHWA,NBasis,NBasis);
        	  FockHWB = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outFockHWB,NBasis,NBasis);
        	  GibbsYA = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYA,NBasis,NBasis);
        	  GibbsYB = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYB,NBasis,NBasis);

        	  GibbsYKernel1A = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYKernel1A,NBasis,NBasis);
        	  GibbsYKernel1B = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYKernel1B,NBasis,NBasis);
        	  GibbsYKernel2A = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYKernel2A,NBasis,NBasis);
        	  GibbsYKernel2B = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYKernel2B,NBasis,NBasis);
        	  /*
        	  GibbsYKernel1A = Eigen::Map<Eigen::Matrix<std::complex< double >,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYKernel1A,NBasis,NBasis);
        	  GibbsYKernel1B = Eigen::Map<Eigen::Matrix<std::complex< double >,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYKernel1B,NBasis,NBasis);
        	  GibbsYKernel2A = Eigen::Map<Eigen::Matrix<std::complex< double >,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYKernel2A,NBasis,NBasis);
        	  GibbsYKernel2B = Eigen::Map<Eigen::Matrix<std::complex< double >,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >(outGibbsYKernel2B,NBasis,NBasis);
        	  */

        	  cout << "START printEigenmatrixCleanFmt(FockHWA)\n";
        	  printEigenmatrixCleanFmt(FockHWA);
        	  cout << "END printEigenmatrixCleanFmt(FockHWA)\n";
        	  cout << "START printEigenmatrixCleanFmt(FockHWB)\n";
        	  printEigenmatrixCleanFmt(FockHWB);
        	  cout << "END printEigenmatrixCleanFmt(FockHWB)\n";
        	  cout << "START printEigenmatrixCleanFmt(GibbsYA)\n";
        	  printEigenmatrixCleanFmt(GibbsYA);
        	  cout << "END printEigenmatrixCleanFmt(GibbsYA)\n";
        	  cout << "START printEigenmatrixCleanFmt(GibbsYB)\n";
        	  printEigenmatrixCleanFmt(GibbsYB);
        	  cout << "END printEigenmatrixCleanFmt(GibbsYB)\n";
        	  cout << "---------------------------------------------------------\n";
        	  cout << "------------------ END LAST ANTC CYCLE ------------------\n";
        	  cout << "---------------------------------------------------------\n";
          }


        //} while (((ediff_rel > conv) || (rms_error > conv)) && (iter < maxiter));
        } while (((ediff_rel > conv) || (Falpha_err > conv) || (iter < 20) ) && (iter < maxiter));
        cout << "EXIT do_ant IF\n";
    }
    cout << "*** =========================== ***\n";
    cout << "***        END ANT loop         ***\n";
    cout << "*** =========================== ***\n";
    cout << "*** =========================== ***\n";
    cout << "***        END SCF loop         ***\n";
    cout << "*** =========================== ***\n";




    printf("\n\n\n");
    printf("********************************************************************************\n");
    printf("** Hartree-Fock energy = %20.12f\n", ehf + enuc);
    printf("********************************************************************************\n");
    printf(" SCF Done:  E(UHF) = %20.12f     A.U. after   %d cycles\n", ehf + enuc,inputNCycles);
    printf("********************************************************************************\n");
    printf("\n\n\n");

    auto Mu = compute_1body_ints<Operator::emultipole2>(obs);
    std::array<double, 3> mu;
    std::array<double, 6> qu;
    for (int xyz = 0; xyz != 3; ++xyz)
      mu[xyz] = -2 *
                D.cwiseProduct(Mu[xyz + 1])
                    .sum();  // 2 = alpha + beta, -1 = electron charge
    for (int k = 0; k != 6; ++k)
      qu[k] = -2 *
              D.cwiseProduct(Mu[k + 4])
                  .sum();  // 2 = alpha + beta, -1 = electron charge
    cout << "** edipole = ";
    std::copy(mu.begin(), mu.end(),
              std::ostream_iterator<double>(cout, " "));
    cout << std::endl;
    cout << "** equadrupole = ";
    std::copy(qu.begin(), qu.end(),
              std::ostream_iterator<double>(cout, " "));
    cout << std::endl;



    cout << "---------------------------------------------------------\n";
    cout << "----------------- AFTER LAST ANTC CYCLE -----------------\n";
    cout << "--------------------- COMPUTE FORCE ---------------------\n";
    cout << "---------------------------------------------------------\n";

    {  // compute force
#if LIBINT2_DERIV_ONEBODY_ORDER
      // compute 1-e forces
      Matrix F1 = Matrix::Zero(atoms.size(), 3);
      Matrix F_Pulay = Matrix::Zero(atoms.size(), 3);
      cout << "BEFORE one-body contributions to the forces\n";
      //////////
      // one-body contributions to the forces
      //////////
      auto T1 = compute_1body_ints_deriv<Operator::kinetic>(1, obs, atoms);
      auto V1 = compute_1body_ints_deriv<Operator::nuclear>(1, obs, atoms);
      auto H01 = T1; // WILL CONTAIN 1-E HAMILTONIAN, SO I HAVE TO ADD VECTOR OF MATRICES V1.
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          H01[i] += V1[i];
          // NUMBER 2 SEEMS TO BE DUE TO SPIN AND MAY NEED TO BE REMOVED.
          auto force = 2 * (T1[i] + V1[i]).cwiseProduct(D).sum();
          F1(atom, xyz) += force;
        }
      }


      // PULAY FORCES ARE RELATED TO OVERLAP MATRIX ELEMENTS AND ITS DERIVATIVES.
      cout << "PULAY FORCES ARE RELATED TO OVERLAP MATRIX ELEMENTS AND ITS DERIVATIVES.\n";
      cout << "BEFORE Pulay force\n";
      //////////
      // Pulay force
      //////////
      //---------------------------------
      // evals seems to contain the Kohn-Sham Eigenvalues.
      // evals_occ seems to contain the ndocc occuppied Kohn-Sham Eigenvalues.
      // INSTEAD OF SUMMING THE ndocc LOWEST K-S EIGENVALUES,
      // THE INTEGRAL OF THE FOCK OPERATOR SEEMS TO BE NECESSARY.
      // THE FOCK OPERATOR IS REAL, HAS REAL EIGENVALUES.
      // USING THE EIGENVALUES SEEMS TO BE STILL LEGIT. NOOO!!!
      // THE INTEGRATION OF THE FOCK OPERATOR WITH A FROZEN DENSITY NEEDS TO BE PERFORMED.
      // THE INTEGRAL OF THE FOCK OPERATOR REPLACES THE MATRIX W BELOW.
      // 5 LINES BELOW EXTRACTED FROM DENSITY MATRIX CALCULATION.
      //Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(X.transpose() * F_diis * X);
      //evals = eig_solver.eigenvalues();
      //auto C = X * eig_solver.eigenvectors();
      //C_occ = C.leftCols(ndocc);
      //D = C_occ * C_occ.transpose();
      //---------------------------------
      cout << "BEFORE evals_occ\n";
      // orbital energy density
      //DiagonalMatrix evals_occ(evals.topRows(ndocc)); // ORIGINAL, COMMENTED.
      //DiagonalMatrix evals_occ(0.5*evals.topRows(naocc) + 0.5*evals.topRows(nbocc));
      cout << "BEFORE W\n";
      //Matrix W = C_occ * evals_occ * C_occ.transpose(); // ORIGINAL, COMMENTED.
      Matrix W = FockHWA + FockHWB;
      cout << "START printEigenmatrixCleanFmt(W)\n";
      printEigenmatrixCleanFmt(W);
      cout << "END printEigenmatrixCleanFmt(W)\n";
      cout << "BEFORE S1\n";
      auto S1 = compute_1body_ints_deriv<Operator::overlap>(1, obs, atoms);
      cout << "BEFORE F_Pulay\n";
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          //auto force = 2 * S1[i].cwiseProduct(W).sum(); // ORIGINAL, COMMENTED. 2 SEEMS TO BE DUE TO SPIN.
          auto force = S1[i].cwiseProduct(W).sum();
          F_Pulay(atom, xyz) -= force;
        }
      }

      cout << "** 1-body forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << F1(atom, xyz) << " ";
      cout << std::endl;

      cout << "START printEigenmatrixCleanFmt(F1)\n";
      printEigenmatrixCleanFmt(F1);
      cout << "END printEigenmatrixCleanFmt(F1)\n";

      cout << "** Pulay forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz)
          cout << F_Pulay(atom, xyz) << " ";
      cout << std::endl;

      cout << "START printEigenmatrixCleanFmt(F_Pulay)\n";
      printEigenmatrixCleanFmt(F_Pulay);
      cout << "END printEigenmatrixCleanFmt(F_Pulay)\n";

#endif  // LIBINT2_DERIV_ONEBODY_ORDER







std::vector<std::vector<Matrix> > G1;
auto natoms = atoms.size();
G1.resize(3*natoms); // natoms, atoms.size()
Matrix F2 = Matrix::Zero(atoms.size(), 3);
for(size_t i=0;i<3*natoms;i++) {
	G1[i].resize(2);
	G1[i][0].resize(Calpha.rows(),Calpha.cols());
	G1[i][1].resize(Calpha.rows(),Calpha.cols());
}


if (not do_density_fitting) {

#if LIBINT2_DERIV_ERI_ORDER
      cout << "BEFORE compute 2-e forces\n";
      // compute 2-e forces
      //Matrix F2 = Matrix::Zero(atoms.size(), 3);

      //////////
      // two-body contributions to the forces
      //////////
      //auto G1 = compute_2body_fock_deriv<1>(obs, atoms, D);
      //auto G1 = compute_2body_fock_uhf_deriv<1>(obs, atoms, Dalpha, Dbeta);
      G1 = compute_2body_fock_uhf_deriv<1>(obs, atoms, Dalpha, Dbeta);
      cout << " EXITED compute_2body_fock_uhf_deriv<1>(obs, atoms, Dalpha, Dbeta); ";
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
          //auto force = G1[i].cwiseProduct(D).sum();
          //cout << "START printEigenmatrixCleanFmt(G1["<<i<<"][0])\n";
          //printEigenmatrixCleanFmt(G1[i][0]);
          //cout << "END printEigenmatrixCleanFmt(G1["<<i<<"][0])\n";
          //cout << "START printEigenmatrixCleanFmt(G1["<<i<<"][1])\n";
          //printEigenmatrixCleanFmt(G1[i][1]);
          //cout << "END printEigenmatrixCleanFmt(G1["<<i<<"][1])\n";
          auto force = G1[i][0].cwiseProduct(Dalpha).sum();
          force += G1[i][1].cwiseProduct(Dbeta).sum();
          F2(atom, xyz) += force;
        }
      }

      cout << "** 2-body forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << F2(atom, xyz) << " ";
      cout << std::endl;

      cout << "START printEigenmatrixCleanFmt(F2)\n";
      printEigenmatrixCleanFmt(F2);
      cout << "END printEigenmatrixCleanFmt(F2)\n";

#endif

}







/*
#if LIBINT2_DERIV_ERI_ORDER
      cout << "BEFORE compute 2-e forces\n";
      // compute 2-e forces
      //Matrix F2 = Matrix::Zero(atoms.size(), 3);
      Matrix FJ = Matrix::Zero(atoms.size(), 3);

      //////////
      // two-body contributions to the forces
      //////////
      //auto G1 = compute_2body_fock_deriv<1>(obs, atoms, D);
      //auto G1 = compute_2body_JK_uhf_deriv<1>(obs, atoms, Dalpha, Dbeta);
      auto JK1_add = compute_2body_JK_uhf_deriv<1>(obs, Dalpha_diff, Dbeta_diff, precision_F, K);
      if (debugflag){
        cout << "START printEigenmatrix(JK_add)\n";
        //printEigenmatrix(JK1_add);
        cout << "END printEigenmatrix(JK_add)\n";
      }

      auto J1 = JK1_add[0];
      auto K1alpha = JK1_add[1];
      auto K1beta = JK1_add[2];
      cout << " EXITED compute_2body_JK_uhf_deriv<1>(obs, atoms, Dalpha, Dbeta); ";
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
          //auto force = G1[i].cwiseProduct(D).sum();
          //cout << "START printEigenmatrixCleanFmt(G1["<<i<<"][0])\n";
          //printEigenmatrixCleanFmt(G1[i][0]);
          //cout << "END printEigenmatrixCleanFmt(G1["<<i<<"][0])\n";
          //cout << "START printEigenmatrixCleanFmt(G1["<<i<<"][1])\n";
          //printEigenmatrixCleanFmt(G1[i][1]);
          //cout << "END printEigenmatrixCleanFmt(G1["<<i<<"][1])\n";
          //auto force = G1[i][0].cwiseProduct(Dalpha).sum();
          auto forceJ = J1[i].cwiseProduct(Dalpha + Dbeta).sum();
          auto forceK = K1alpha[i].cwiseProduct(Dalpha).sum();
          forceK += K1beta[i][2].cwiseProduct(Dbeta).sum();
          //F2(atom, xyz) += force;
          FJ(atom, xyz) += forceJ;
        }
      }

      cout << "** 2-body forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << FJ(atom, xyz) << " ";
      cout << std::endl;

      cout << "START printEigenmatrixCleanFmt(FJ)\n";
      printEigenmatrixCleanFmt(FJ);
      cout << "END printEigenmatrixCleanFmt(FJ)\n";

#endif
*/



   #if HAVE_DENSITY_FITTING
         else {  // do DF
          if (debugflag)cout << "if (do_density_fitting)\n";
       	  if (debugflag)cout << "INSIDE HAVE_DENSITY_FITTING\n";
       	  if (debugflag)cout << "WHILE CALCULATING DFT FORCES\n";
       	  std::vector<>
          //   //F = H + dffockengine->compute_2body_fock_dfC(C_occ);
       	  //Falpha = H + dffockengine->compute_2body_fock_dfC(Calpha_occ);
       	  //Fbeta = H + dffockengine->compute_2body_fock_dfC(Cbeta_occ);
       	  //F = {Falpha,Fbeta};
       	  //cout<<"Press ENTER to continue...";cin.get();
         }
   #else
         else {  // yes do DF
        	 if (debugflag)cout << "if (do_density_fitting)\n";
        	 if (debugflag)cout << "INSIDE NOT HAVE_DENSITY_FITTING\n";
        	 if (debugflag)cout << "WHILE CALCULATING DFT FORCES\n";

#if LIBINT2_DERIV_ERI_ORDER
      cout << "BEFORE compute 2-e forces\n";
      // compute 2-e forces
      //Matrix F2 = Matrix::Zero(atoms.size(), 3);

      //////////
      // two-body contributions to the forces
      //////////
      // VARIABLES QUE ES NECESARIO ACTUALIZAR.
      cout<<"VARIABLES QUE ES NECESARIO ACTUALIZAR.\n";
      const auto precision_F = std::min(
              			  std::min(1e-3 / XtX_condition_number, 1e-7),
      					  std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon()));

      if (debugflag){
       	cout << "rms_error = " << rms_error<<"\n";
        cout << "1e-3 / XtX_condition_number = " << 1e-3 / XtX_condition_number<<"\n";
        cout << "XtX_condition_number = " << XtX_condition_number<<"\n";
        cout << "std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon()) = " << std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon())<<"\n";
        cout << "INSIDE if (not do_density_fitting) BEFORE compute_2body_fock\n";
        cout << "START printEigenmatrix(D_diff)\n";
        //printEigenmatrix(D_diff);
        cout << "END printEigenmatrix(D_diff)\n";
        cout << "START printEigenmatrix(K)\n";
        //printEigenmatrix(K);
        cout << "END printEigenmatrix(K)\n";
        cout << "START printEigenmatrix(F)\n";
        //printEigenmatrix(F);
        cout << "END printEigenmatrix(F)\n";
        cout << "precision_F = " << precision_F << "\n";
      }

      Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_alpha(X.transpose() * Falpha * X);
      Eigen::SelfAdjointEigenSolver<Matrix> eig_solver_beta(X.transpose() * Fbeta * X);
      //auto Calpha = X * eig_solver_alpha.eigenvectors();
      //auto Cbeta = X * eig_solver_beta.eigenvectors();
      Matrix Calpha = X * eig_solver_alpha.eigenvectors();
      Matrix Cbeta = X * eig_solver_beta.eigenvectors();
      evalsalpha = eig_solver_alpha.eigenvalues();
      evalsbeta = eig_solver_beta.eigenvalues();

      // compute density, D = C(occ).C(occ)T
      Calpha_occ = Calpha.leftCols(naocc);
      //Dalpha = Calpha_occ * Calpha_occ.transpose();
      Cbeta_occ = Cbeta.leftCols(nbocc);
      //Dbeta = Cbeta_occ * Cbeta_occ.transpose();
      //D = Dalpha + Dbeta;
      //F += compute_2body_fock(obs, D_diff, precision_F, K);
      auto JK_add = compute_2body_JK_uhf(obs, Dalpha_diff, Dbeta_diff, precision_F, K);
      if (debugflag){
    	  cout << "START printEigenmatrix(JK_add)\n";
    	  //printEigenmatrix(F_add);
    	  cout << "END printEigenmatrix(JK_add)\n";
      }

      auto J = JK_add[0];
      auto Kalpha = JK_add[1];
      auto Kbeta = JK_add[2];

      Matrix XCalpha = Matrix::Zero(NBasis,NBasis);
      Matrix XCbeta = Matrix::Zero(NBasis,NBasis);

      if (debugflag)cout << "CALLING DFT FROM ERKALE\n";

	  cout << "START printEigenmatrix(evalsalpha)\n";
	  printEigenmatrix(evalsalpha);
	  cout << "END printEigenmatrix(evalsalpha)\n";

      //double Ea = 0.0;
      //double Eb = 0.0;
      //std::vector<double> Ea(evalsalpha.diagonal().data(),evalsalpha.diagonal().data() + evalsalpha.diagonal().size());
      //std::vector<double> Eb(evalsbeta.diagonal().data(),evalsbeta.diagonal().data() + evalsbeta.diagonal().size());
	  // evalsalpha y evalsbeta NO SON MATRICES DIAGONALES, SON VECTORES.
	  std::vector<double> Ea(evalsalpha.data(),evalsalpha.data() + evalsalpha.size());
	  std::vector<double> Eb(evalsbeta.data(),evalsbeta.data() + evalsbeta.size());

      cout <<"evalsalpha.size() = "<<evalsalpha.size()<<"\n";

      cout <<"Ea";
      printStdVector(Ea);
      cout <<"Eb";
      printStdVector(Eb);
      cout<<"END VARIABLES QUE ES NECESARIO ACTUALIZAR.\n";
      // END VARIABLES QUE ES NECESARIO ACTUALIZAR.


      //auto G1 = compute_2body_fock_deriv<1>(obs, atoms, D);
      //auto G1 = compute_2body_fockj_uhf_deriv<1>(obs, atoms, Dalpha, Dbeta);
      //wrap_Fock_UDFT(ErkaleSCF,erkalebasis,erkalejobname,NBasis,Ea,Eb,naocc,nbocc,Calpha,Cbeta,Calpha_occ,Cbeta_occ,Falpha,Fbeta,Dalpha,Dbeta,J,Kalpha,Kbeta,XCalpha,XCbeta);
      //auto G1Erk = wrap_force_UDFT(ErkaleSCF,erkalebasis,erkalejobname,NBasis,Ea,Eb,naocc,nbocc,Calpha,Cbeta,Calpha_occ,Cbeta_occ,Halpha,Hbeta,Dalpha,Dbeta,J,Kalpha,Kbeta,XCalpha,XCbeta);
      //wrap_Fock_UDFT(ErkaleSCF,erkalebasis,erkalejobname,NBasis,Ea,Eb,naocc,nbocc,Calpha,Cbeta,Calpha_occ,Cbeta_occ,Falpha,Fbeta,Dalpha,Dbeta,J,Kalpha,Kbeta,XCalpha,XCbeta);
      // PASAMOS Falpha y Fbeta por Halpha y Hbeta igual que en wrap_Fock_UDFT.
      G1 = wrap_force_UDFT(ErkaleSCF,erkalebasis,erkalesettings,erkalejobname,NBasis,Ea,Eb,naocc,nbocc,S,Calpha,Cbeta,Calpha_occ,Cbeta_occ,Falpha,Fbeta,Dalpha,Dbeta,J,Kalpha,Kbeta,XCalpha,XCbeta);
      cout << " EXITED compute_2body_fock_uhf_deriv<1>(obs, atoms, Dalpha, Dbeta); \n";
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
          //auto force = G1[i].cwiseProduct(D).sum();
          //cout << "START printEigenmatrixCleanFmt(G1["<<i<<"][0])\n";
          //printEigenmatrixCleanFmt(G1[i][0]);
          //cout << "END printEigenmatrixCleanFmt(G1["<<i<<"][0])\n";
          //cout << "START printEigenmatrixCleanFmt(G1["<<i<<"][1])\n";
          //printEigenmatrixCleanFmt(G1[i][1]);
          //cout << "END printEigenmatrixCleanFmt(G1["<<i<<"][1])\n";
          //auto force = G1[i][0].cwiseProduct(Dalpha).sum();
          //force += G1[i][1].cwiseProduct(Dbeta).sum();
          auto force = G1[i][0].trace();
          force += G1[i][1].trace();
          //auto force = G1Erk[i]; // VER QUE EL INCREMENTAL i VA DESDE 1 HASTA 3*basis.get_Nnuc();
          //auto force = G1[i]; // VER QUE EL INCREMENTAL i VA DESDE 1 HASTA 3*basis.get_Nnuc();
          F2(atom, xyz) += force;
        }
      }

      cout << "\nffull EN scf-force DEBERÍA COINCIDIR CON F2 en HartreeFockClass 2-body Erkale forces\n";
      cout << "** 2-body Erkale forces = \n";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << F2(atom, xyz) << " ";
      cout << std::endl;

      cout << "START printEigenmatrixCleanFmt(F2)\n";
      printEigenmatrixCleanFmt(F2);
      cout << "END printEigenmatrixCleanFmt(F2)\n";






/*
      cout << "BEFORE compute nuclear repulsion forces\n";
      // compute nuclear repulsion forces
      Matrix FN = Matrix::Zero(atoms.size(), 3);
      //////////
      // nuclear repulsion contribution to the forces
      //////////

      cout << "** nuclear repulsion forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << FN(atom, xyz) << " ";
      cout << std::endl;

      cout << "BEFORE Hartree-Fock forces\n";
      //auto F = F1 + F_Pulay + F2 + FN + F_Gibbs;
      //Matrix F = Matrix::Zero(atoms.size(), 3);
      Matrix F = F1 + F_Pulay + F2 + FN + F_Gibbs;
      cout << "** Hartree-Fock forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << F(atom, xyz) << " ";
      cout << std::endl;

      cout << "START printEigenmatrixCleanFmt(F)\n";
      printEigenmatrixCleanFmt(F);
      cout << "END printEigenmatrixCleanFmt(F)\n";
      cout << "START writeEigenMatrixToCSVfile(string name, Matrix F)\n";
      //char filename[strlen(basename)+strlen(basisname)+4];
      char outforcesfile[strlen(fullname)+7];
      strcpy(outforcesfile,fullname);
      strcat(outforcesfile,".forces");
      writeEigenMatrixToCSVfile(outforcesfile,F);
      cout << "END writeEigenMatrixToCSVfile(string name, Matrix F)\n";
*/



#endif // END #if LIBINT2_DERIV_ERI_ORDER


         }

   #endif // END #if HAVE_DENSITY_FITTING







// if support 1-e and 2-e derivatives compute the derivatives of the Gibbs Free Energy.
#if LIBINT2_DERIV_ONEBODY_ORDER && LIBINT2_DERIV_ERI_ORDER
      cout << "BEFORE compute Gibbs Free Energy forces\n";
      Matrix F_Gibbs = Matrix::Zero(atoms.size(), 3);

      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {

          /*
          cout << "START printEigenmatrixCleanFmt(H01["<<i<<"])\n";
          printEigenmatrixCleanFmt(H01[i]);
          cout << "END printEigenmatrixCleanFmt(H01["<<i<<"])\n";
          */

          cout << "START printEigenmatrixCleanFmt(G1["<<i<<"][0])\n";
          printEigenmatrixCleanFmt(G1[i][0]);
          cout << "END printEigenmatrixCleanFmt(G1["<<i<<"][0])\n";

          /*
          cout << "START printEigenmatrixCleanFmt(G1["<<i<<"][1])\n";
          printEigenmatrixCleanFmt(G1[i][1]);
          cout << "END printEigenmatrixCleanFmt(G1["<<i<<"][1])\n";

          cout << "START printEigenmatrixCleanFmt(S1["<<i<<"])\n";
          printEigenmatrixCleanFmt(S1[i]);
          cout << "END printEigenmatrixCleanFmt(S1["<<i<<"])\n";
          */

          Matrix S1piv = S1[i];
          Matrix Halphapiv = H01[i] + G1[i][0];
          Matrix Hbetapiv = H01[i] + G1[i][1];

          cout << "i="<<i<<"\n";
          //cout << "START printEigenmatrixCleanFmt(Halphapiv);\n";
          //printEigenmatrixCleanFmt(Halphapiv);
          //cout << "END printEigenmatrixCleanFmt(Halphapiv);\n";

          // EN LOS PRODUCTOS DE MATRICES NO PODEMOS DEJAR ESPACIO. SI DEJAMOS ESPACIO HACE cwiseProduct Y NO PRODUCTO ORDINARIO.
          //auto force = 2 * S1[i].cwiseProduct(W).sum(); // ORIGINAL, COMMENTED. 2 SEEMS TO BE DUE TO SPIN.
          auto force = ((GibbsYKernel1A*S.inverse()*S1piv*S.inverse()).cwiseProduct(Dalpha)).sum();
          force += ((GibbsYKernel1B*S.inverse()*S1piv*S.inverse()).cwiseProduct(Dbeta)).sum();
          force += ((GibbsYKernel2A*S.inverse()*Halphapiv*S.inverse()).cwiseProduct(Dalpha)).sum();
          force += ((GibbsYKernel2B*S.inverse()*Hbetapiv*S.inverse()).cwiseProduct(Dbeta)).sum();
          F_Gibbs(atom, xyz) -= force;
        }
      }

      cout << "START printEigenmatrixCleanFmt(S)\n";
      printEigenmatrixCleanFmt(S);
      cout << "END printEigenmatrixCleanFmt(S)\n";

      cout << "START printEigenmatrixCleanFmt(GibbsYKernel1A)\n";
      printEigenmatrixCleanFmt(GibbsYKernel1A);
      cout << "END printEigenmatrixCleanFmt(GibbsYKernel1A)\n";

      cout << "START printEigenmatrixCleanFmt(GibbsYKernel1B)\n";
      printEigenmatrixCleanFmt(GibbsYKernel1B);
      cout << "END printEigenmatrixCleanFmt(GibbsYKernel1B)\n";

      cout << "START printEigenmatrixCleanFmt(GibbsYKernel2A)\n";
      printEigenmatrixCleanFmt(GibbsYKernel2A);
      cout << "END printEigenmatrixCleanFmt(GibbsYKernel2A)\n";

      cout << "START printEigenmatrixCleanFmt(GibbsYKernel2B)\n";
      printEigenmatrixCleanFmt(GibbsYKernel2B);
      cout << "END printEigenmatrixCleanFmt(GibbsYKernel2B)\n";

      cout << "START printEigenmatrixCleanFmt(Dalpha)\n";
      printEigenmatrixCleanFmt(Dalpha);
      cout << "END printEigenmatrixCleanFmt(Dalpha)\n";

      cout << "START printEigenmatrixCleanFmt(Dbeta)\n";
      printEigenmatrixCleanFmt(Dbeta);
      cout << "END printEigenmatrixCleanFmt(Dbeta)\n";

      cout << "\n\n";
      cout << "** 2-body forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << F_Gibbs(atom, xyz) << " ";
      cout << std::endl;

      cout << "\n\n";
      cout << "********************************************************\n";
      cout << "********************************************************\n";
      cout << "********************************************************\n";
      cout << "START printEigenmatrixCleanFmt(F_Gibbs)\n";
      printEigenmatrixCleanFmt(F_Gibbs);
      cout << "END printEigenmatrixCleanFmt(F_Gibbs)\n";
      cout << "********************************************************\n";
      cout << "********************************************************\n";
      cout << "********************************************************\n";
      cout << "\n\n";

      cout << "AFTER compute Gibbs Free Energy forces\n";
      cout << "\n\n";


      cout << "\n\n";
      cout << "********************************************************\n";
      cout << "********************************************************\n";
      cout << "********************************************************\n";
      cout << "START FOCK MATRIX printEigenmatrixCleanFmt(Falpha)\n";
      printEigenmatrixCleanFmt(Falpha);
      cout << "END FOCK MATRIX printEigenmatrixCleanFmt(Falpha)\n";
      cout << "********************************************************\n";
      cout << "********************************************************\n";
      cout << "********************************************************\n";
      cout << "\n\n";
#endif

// if support 1-e and 2-e derivatives compute nuclear repulsion force and the
// total force
#if LIBINT2_DERIV_ONEBODY_ORDER && LIBINT2_DERIV_ERI_ORDER
      cout << "BEFORE compute nuclear repulsion forces\n";
      // compute nuclear repulsion forces
      Matrix FN = Matrix::Zero(atoms.size(), 3);
      //////////
      // nuclear repulsion contribution to the forces
      //////////
      for (auto a1 = 1; a1 != atoms.size(); ++a1) {
        const auto& atom1 = atoms[a1];
        for (auto a2 = 0; a2 < a1; ++a2) {
          const auto& atom2 = atoms[a2];

          auto x12 = atom1.x - atom2.x;
          auto y12 = atom1.y - atom2.y;
          auto z12 = atom1.z - atom2.z;
          auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
          auto r12 = sqrt(r12_2);
          auto r12_3 = r12 * r12_2;

          auto z1z2_over_r12_3 =
              atom1.atomic_number * atom2.atomic_number / r12_3;

          auto fx = -x12 * z1z2_over_r12_3;
          auto fy = -y12 * z1z2_over_r12_3;
          auto fz = -z12 * z1z2_over_r12_3;
          FN(a1, 0) += fx;
          FN(a1, 1) += fy;
          FN(a1, 2) += fz;
          FN(a2, 0) -= fx;
          FN(a2, 1) -= fy;
          FN(a2, 2) -= fz;
        }
      }

      cout << "** nuclear repulsion forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << FN(atom, xyz) << " ";
      cout << std::endl;

      cout << "BEFORE Hartree-Fock forces\n";
      //auto F = F1 + F_Pulay + F2 + FN + F_Gibbs;
      //Matrix F = Matrix::Zero(atoms.size(), 3);
      Matrix F = F1 + F_Pulay + F2 + FN + F_Gibbs;
      cout << "** Hartree-Fock forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) cout << F(atom, xyz) << " ";
      cout << std::endl;

      cout << "START printEigenmatrixCleanFmt(F)\n";
      printEigenmatrixCleanFmt(F);
      cout << "END printEigenmatrixCleanFmt(F)\n";
      cout << "START writeEigenMatrixToCSVfile(string name, Matrix F)\n";
      //char filename[strlen(basename)+strlen(basisname)+4];
      char outforcesfile[strlen(fullname)+7];
      strcpy(outforcesfile,fullname);
      strcat(outforcesfile,".forces");
      writeEigenMatrixToCSVfile(outforcesfile,F);
      cout << "END writeEigenMatrixToCSVfile(string name, Matrix F)\n";


#endif
    }

    {  // compute hessian
      cout << "BEFORE compute hessian\n";
      const auto ncoords = 3 * atoms.size();
      // # of elems in upper triangle
      const auto nelem =  ncoords * (ncoords+1) / 2;
#if LIBINT2_DERIV_ONEBODY_ORDER > 1
      // compute 1-e hessian
      Matrix H1 = Matrix::Zero(ncoords, ncoords);
      Matrix H_Pulay = Matrix::Zero(ncoords, ncoords);
      cout << "BEFORE one-body contributions to the hessian\n";
      //////////
      // one-body contributions to the hessian
      //////////
      auto T2 = compute_1body_ints_deriv<Operator::kinetic>(2, obs, atoms);
      auto V2 = compute_1body_ints_deriv<Operator::nuclear>(2, obs, atoms);
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          auto hess = 2 * (T2[i] + V2[i]).cwiseProduct(D).sum();
          H1(row, col) += hess;
        }
      }

      cout << "BEFORE Pulay hessian\n";
      //////////
      // Pulay hessian
      //////////
      // orbital energy density
      //DiagonalMatrix evals_occalpha(evals.topRows(naocc)); // ORIGINAL CODE. COMMENTED.
      //DiagonalMatrix evals_occbeta(evals.topRows(nbocc)); // ORIGINAL CODE. COMMENTED.
      //Matrix W = Calpha_occ * evals_occalpha * Calpha_occ.transpose() + Cbeta_occ * evals_occbeta * Cbeta_occ.transpose();
      Matrix W = FockHWA + FockHWB;
      cout << "START printEigenmatrixCleanFmt(W)\n";
      printEigenmatrixCleanFmt(W);
      cout << "END printEigenmatrixCleanFmt(W)\n";
      auto S2 = compute_1body_ints_deriv<Operator::overlap>(2, obs, atoms);
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          //auto hess = 2 * S2[i].cwiseProduct(W).sum(); // ORIGINAL CODE. COMMENTED. 2 SEEMS TO BE DUE TO SPIN.
          auto hess = S2[i].cwiseProduct(W).sum();
          H_Pulay(row, col) -= hess;
        }
      }

      cout << "** 1-body hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << H1(row, col) << " ";
        }
      }
      cout << std::endl;

      cout << "** Pulay hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << H_Pulay(row, col) << " ";
        }
      }
      cout << std::endl;
#endif  // LIBINT2_DERIV_ONEBODY_ORDER > 1

#if LIBINT2_DERIV_ERI_ORDER > 1
      cout << "BEFORE compute 2-e hessian\n";
      // compute 2-e forces
      Matrix H2 = Matrix::Zero(ncoords, ncoords);

      //////////
      // two-body contributions to the forces
      //////////
      //auto G2 = compute_2body_fock_deriv<2>(obs, atoms, D);
      auto G2 = compute_2body_fock_uhf_deriv<2>(obs, atoms, Dalpha, Dbeta);
      cout << " EXITED compute_2body_fock_uhf_deriv<2>(obs, atoms, Dalpha, Dbeta); ";

      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
          //cout << "START printEigenmatrixCleanFmt(G2["<<i<<"][0])\n";
          //printEigenmatrixCleanFmt(G2[i][0]);
          //cout << "END printEigenmatrixCleanFmt(G2["<<i<<"][0])\n";
          //cout << "START printEigenmatrixCleanFmt(G2["<<i<<"][1])\n";
          //printEigenmatrixCleanFmt(G2[i][1]);
          //cout << "END printEigenmatrixCleanFmt(G2["<<i<<"][1])\n";
          auto hess = G2[i][0].cwiseProduct(Dalpha).sum();
          hess += G2[i][1].cwiseProduct(Dbeta).sum();
          H2(row, col) += hess;
        }
      }

      cout << "** 2-body hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << H2(row, col) << " ";
        }
      }
      cout << std::endl;
#endif


// if support 1-e and 2-e 2nd derivatives compute the Gibbs Free Energy hessian
//  HESSIAN IS NOT NECESARY FOR BERNY ALGORITHM.
/*
#if LIBINT2_DERIV_ONEBODY_ORDER > 1 && LIBINT2_DERIV_ERI_ORDER > 1
      cout << "BEFORE compute Gibbs Free Energy hessian\n";
      Matrix F_Gibbs = Matrix::Zero(atoms.size(), 3);

      // 1-body first derivative operators.
      auto T1 = compute_1body_ints_deriv<Operator::kinetic>(1, obs, atoms);
      auto V1 = compute_1body_ints_deriv<Operator::nuclear>(1, obs, atoms);

      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          //auto Halphapiv = H01[i] + G1[i][0];
          //auto Hbetapiv = H01[i] + G1[i][1];
          auto Halphapiv = T1[i] + V1[i] + G1[i][0];
          auto Hbetapiv = T1[i] + V1[i] + G1[i][1];
          //auto force = 2 * S1[i].cwiseProduct(W).sum(); // ORIGINAL, COMMENTED. 2 SEEMS TO BE DUE TO SPIN.
          auto force = (GibbsYKernel1A * S.inverse() * S1[i] * S.inverse()).cwiseProduct(Dalpha).sum();
          force += (GibbsYKernel1B * S.inverse() * S1[i] * S.inverse()).cwiseProduct(Dbeta).sum();
          force += (GibbsYKernel2A * S.inverse() * Halphapiv * S.inverse()).cwiseProduct(Dalpha).sum();
          force += (GibbsYKernel2B * S.inverse() * Hbetapiv * S.inverse()).cwiseProduct(Dbeta).sum();
          //F_Gibbs(atom, xyz) += force;
        }
      }
      cout << "START printEigenmatrixCleanFmt(F_Gibbs)\n";
      printEigenmatrixCleanFmt(F_Gibbs);
      cout << "END printEigenmatrixCleanFmt(F_Gibbs)\n";

      cout << "AFTER compute Gibbs Free Energy hessian\n";
#endif
*/


// if support 1-e and 2-e 2nd derivatives compute nuclear repulsion hessian and
// the total hessian
#if LIBINT2_DERIV_ONEBODY_ORDER > 1 && LIBINT2_DERIV_ERI_ORDER > 1
      cout << "BEFORE compute nuclear repulsion hessian\n";
      // compute nuclear repulsion hessian
      // NB only the upper triangle is computed!!!
      Matrix HN = Matrix::Zero(ncoords, ncoords);

      //////////
      // nuclear repulsion contribution to the hessian
      //////////
      for (auto a1 = 1; a1 != atoms.size(); ++a1) {
        const auto& atom1 = atoms[a1];
        for (auto a2 = 0; a2 < a1; ++a2) {
          const auto& atom2 = atoms[a2];

          auto x12 = atom1.x - atom2.x;
          auto y12 = atom1.y - atom2.y;
          auto z12 = atom1.z - atom2.z;
          auto x12_2 = x12 * x12;
          auto y12_2 = y12 * y12;
          auto z12_2 = z12 * z12;
          auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
          auto r12 = sqrt(r12_2);
          auto r12_5 = r12 * r12_2 * r12_2;

          auto z1z2_over_r12_5 =
              atom1.atomic_number * atom2.atomic_number / r12_5;

          HN(3*a1 + 0, 3*a1 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
          HN(3*a1 + 1, 3*a1 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
          HN(3*a1 + 2, 3*a1 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
          HN(3*a1 + 0, 3*a1 + 1) += z1z2_over_r12_5 * (3*x12*y12);
          HN(3*a1 + 0, 3*a1 + 2) += z1z2_over_r12_5 * (3*x12*z12);
          HN(3*a1 + 1, 3*a1 + 2) += z1z2_over_r12_5 * (3*y12*z12);

          HN(3*a2 + 0, 3*a2 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
          HN(3*a2 + 1, 3*a2 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
          HN(3*a2 + 2, 3*a2 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
          HN(3*a2 + 0, 3*a2 + 1) += z1z2_over_r12_5 * (3*x12*y12);
          HN(3*a2 + 0, 3*a2 + 2) += z1z2_over_r12_5 * (3*x12*z12);
          HN(3*a2 + 1, 3*a2 + 2) += z1z2_over_r12_5 * (3*y12*z12);

          HN(3*a2 + 0, 3*a1 + 0) -= z1z2_over_r12_5 * (3*x12_2 - r12_2);
          HN(3*a2 + 1, 3*a1 + 1) -= z1z2_over_r12_5 * (3*y12_2 - r12_2);
          HN(3*a2 + 2, 3*a1 + 2) -= z1z2_over_r12_5 * (3*z12_2 - r12_2);
          HN(3*a2 + 1, 3*a1 + 0) -= z1z2_over_r12_5 * (3*y12*x12);
          HN(3*a2 + 2, 3*a1 + 0) -= z1z2_over_r12_5 * (3*z12*x12);
          HN(3*a2 + 2, 3*a1 + 1) -= z1z2_over_r12_5 * (3*z12*y12);
          HN(3*a2 + 0, 3*a1 + 1) -= z1z2_over_r12_5 * (3*x12*y12);
          HN(3*a2 + 0, 3*a1 + 2) -= z1z2_over_r12_5 * (3*x12*z12);
          HN(3*a2 + 1, 3*a1 + 2) -= z1z2_over_r12_5 * (3*y12*z12);
        }
      }

      cout << "** nuclear repulsion hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << HN(row, col) << " ";
        }
      }
      cout << std::endl;

      auto H = H1 + H_Pulay + H2 + HN;
      cout << "** Hartree-Fock hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          cout << H(row, col) << " ";
        }
      }
      cout << std::endl;
#endif
    }
    cout << "---------------------------------------------------------\n";
    cout << "----------------- AFTER LAST ANTC CYCLE -----------------\n";
    cout << "------------------- END COMPUTE FORCE -------------------\n";
    cout << "---------------------------------------------------------\n";




    cout << "BEFORE libint2::finalize();\n";
    libint2::finalize();  // done with libint

  }  // end of try block; if any exceptions occurred, report them and exit
     // cleanly

  catch (const char* ex) {
    cerr << "caught exception: " << ex << endl;
    //return 1;
  } catch (std::string& ex) {
    cerr << "caught exception: " << ex << endl;
    //return 1;
  } catch (std::exception& ex) {
    cerr << ex.what() << endl;
    //return 1;
  } catch (...) {
    cerr << "caught unknown exception\n";
    //return 1;
  }

  //return 0; // if int main();

  cout << "---------------------------------------------------------\n";
  cout << "------------- END OF negfyHF.exe EXECUTION --------------\n";
  cout << "---------------------------------------------------------\n";
}

//-------------------------------------------------------------------------------------------
//---- REMOVED FROM RUNANT BEFORE LAST CYCLE ANTC CALL, I PUTTED IT THERE, NOT ORIGINAL -----
//-------------------------------------------------------------------------------------------
/*
{  // compute force
#if LIBINT2_DERIV_ONEBODY_ORDER
  // compute 1-e forces

  cout << "---------------------------------------------------------\n";
  cout << "-------------------- LAST ANTC CYCLE --------------------\n";
  cout << "--------------------- COMPUTE FORCE ---------------------\n";
  cout << "---------------------------------------------------------\n";

  Matrix F1 = Matrix::Zero(atoms.size(), 3);
  Matrix F_Pulay = Matrix::Zero(atoms.size(), 3);
  cout << "BEFORE one-body contributions to the forces\n";
  //////////
  // one-body contributions to the forces
  //////////
  auto T1 = compute_1body_ints_deriv<Operator::kinetic>(1, obs, atoms);
  auto V1 = compute_1body_ints_deriv<Operator::nuclear>(1, obs, atoms);

  Matrix H1x = Matrix::Zero(NBasis,NBasis);
  Matrix H1y = Matrix::Zero(NBasis,NBasis);
  Matrix H1z = Matrix::Zero(NBasis,NBasis);

  for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
    for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
      // NUMBER 2 SEEMS TO BE DUE TO SPIN AND MAY NEED TO BE REMOVED.
      auto force = (T1[i] + V1[i]).cwiseProduct(D).sum(); // COULD ALSO USE Dalpha + Dbeta
      F1(atom, xyz) += force;
    }
  }

  //double * DA = Dalpha.data();
  //double * DB = Dbeta.data();
  auto H1dx = H1x.data();

  // PULAY FORCES ARE RELATED TO OVERLAP MATRIX ELEMENTS AND ITS DERIVATIVES.
  cout << "PULAY FORCES ARE RELATED TO OVERLAP MATRIX ELEMENTS AND ITS DERIVATIVES.\n";
  cout << "BEFORE Pulay force\n";
  //////////
  // Pulay force
  //////////
  //---------------------------------
  // evals seems to contain the Kohn-Sham Eigenvalues.
  // evals_occ seems to contain the ndocc occuppied Kohn-Sham Eigenvalues.
  // INSTEAD OF SUMMING THE ndocc LOWEST K-S EIGENVALUES,
  // THE INTEGRAL OF THE FOCK OPERATOR SEEMS TO BE NECESSARY.
  // THE FOCK OPERATOR IS REAL, HAS REAL EIGENVALUES.
  // USING THE EIGENVALUES SEEMS TO BE STILL LEGIT. NOOO!!!
  // THE INTEGRATION OF THE FOCK OPERATOR WITH A FROZEN DENSITY NEEDS TO BE PERFORMED.
  // THE INTEGRAL OF THE FOCK OPERATOR REPLACES THE MATRIX W BELOW.
  // 5 LINES BELOW EXTRACTED FROM DENSITY MATRIX CALCULATION.
  //Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(X.transpose() * F_diis * X);
  //evals = eig_solver.eigenvalues();
  //auto C = X * eig_solver.eigenvectors();
  //C_occ = C.leftCols(ndocc);
  //D = C_occ * C_occ.transpose();
  //---------------------------------
  cout << "BEFORE evals_occ\n";
  // orbital energy density
  //DiagonalMatrix evals_occ(evals.topRows(ndocc)); // ORIGINAL, COMMENTED.
  //DiagonalMatrix evals_occ(0.5*evals.topRows(naocc) + 0.5*evals.topRows(nbocc));
  cout << "BEFORE W\n";
  //Matrix W = C_occ * evals_occ * C_occ.transpose(); // ORIGINAL, COMMENTED.
  Matrix W = FockHWA + FockHWB;
  cout << "START printEigenmatrixCleanFmt(W)\n";
  printEigenmatrixCleanFmt(W);
  cout << "END printEigenmatrixCleanFmt(W)\n";
  cout << "BEFORE S1\n";
  auto S1 = compute_1body_ints_deriv<Operator::overlap>(1, obs, atoms);
  cout << "BEFORE F_Pulay\n";
  for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
    for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
      //auto force = 2 * S1[i].cwiseProduct(W).sum(); // ORIGINAL, COMMENTED. 2 SEEMS TO BE DUE TO SPIN.
      auto force = S1[i].cwiseProduct(W).sum();
      F_Pulay(atom, xyz) -= force;
    }
  }

  cout << "** 1-body forces = ";
  for (int atom = 0; atom != atoms.size(); ++atom)
    for (int xyz = 0; xyz != 3; ++xyz) cout << F1(atom, xyz) << " ";
  cout << std::endl;
  cout << "** Pulay forces = ";
  for (int atom = 0; atom != atoms.size(); ++atom)
    for (int xyz = 0; xyz != 3; ++xyz)
      cout << F_Pulay(atom, xyz) << " ";
  cout << std::endl;
#endif  // LIBINT2_DERIV_ONEBODY_ORDER

#if LIBINT2_DERIV_ERI_ORDER
  cout << "BEFORE compute 2-e forces\n";
  // compute 2-e forces
  Matrix F2 = Matrix::Zero(atoms.size(), 3);

  //////////
  // two-body contributions to the forces
  //////////
  auto G1 = compute_2body_fock_deriv<1>(obs, atoms, D);
  for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
    for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
      // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
      auto force = G1[i].cwiseProduct(D).sum();
      F2(atom, xyz) += force;
    }
  }

  cout << "** 2-body forces = ";
  for (int atom = 0; atom != atoms.size(); ++atom)
    for (int xyz = 0; xyz != 3; ++xyz) cout << F2(atom, xyz) << " ";
  cout << std::endl;
#endif

// if support 1-e and 2-e derivatives compute nuclear repulsion force and the
// total force
#if LIBINT2_DERIV_ONEBODY_ORDER && LIBINT2_DERIV_ERI_ORDER
  cout << "BEFORE compute nuclear repulsion forces\n";
  // compute nuclear repulsion forces
  Matrix FN = Matrix::Zero(atoms.size(), 3);
  //////////
  // nuclear repulsion contribution to the forces
  //////////
  for (auto a1 = 1; a1 != atoms.size(); ++a1) {
    const auto& atom1 = atoms[a1];
    for (auto a2 = 0; a2 < a1; ++a2) {
      const auto& atom2 = atoms[a2];

      auto x12 = atom1.x - atom2.x;
      auto y12 = atom1.y - atom2.y;
      auto z12 = atom1.z - atom2.z;
      auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
      auto r12 = sqrt(r12_2);
      auto r12_3 = r12 * r12_2;

      auto z1z2_over_r12_3 =
          atom1.atomic_number * atom2.atomic_number / r12_3;

      auto fx = -x12 * z1z2_over_r12_3;
      auto fy = -y12 * z1z2_over_r12_3;
      auto fz = -z12 * z1z2_over_r12_3;
      FN(a1, 0) += fx;
      FN(a1, 1) += fy;
      FN(a1, 2) += fz;
      FN(a2, 0) -= fx;
      FN(a2, 1) -= fy;
      FN(a2, 2) -= fz;
    }
  }

  cout << "** nuclear repulsion forces = ";
  for (int atom = 0; atom != atoms.size(); ++atom)
    for (int xyz = 0; xyz != 3; ++xyz) cout << FN(atom, xyz) << " ";
  cout << std::endl;

  cout << "BEFORE Hartree-Fock forces\n";
  auto F = F1 + F_Pulay + F2 + FN;
  cout << "** Hartree-Fock forces = ";
  for (int atom = 0; atom != atoms.size(); ++atom)
    for (int xyz = 0; xyz != 3; ++xyz) cout << F(atom, xyz) << " ";
  cout << std::endl;

	  cout << "---------------------------------------------------------\n";
	  cout << "-------------------- LAST ANTC CYCLE --------------------\n";
	  cout << "------------------- END COMPUTE FORCE -------------------\n";
	  cout << "---------------------------------------------------------\n";

#endif
}





{  // compute hessian
  cout << "BEFORE compute hessian\n";
  const auto ncoords = 3 * atoms.size();
  // # of elems in upper triangle
  const auto nelem =  ncoords * (ncoords+1) / 2;
#if LIBINT2_DERIV_ONEBODY_ORDER > 1
  // compute 1-e hessian

  cout << "---------------------------------------------------------\n";
  cout << "-------------------- LAST ANTC CYCLE --------------------\n";
  cout << "-------------------- COMPUTE HESSIAN --------------------\n";
  cout << "---------------------------------------------------------\n";

  Matrix H1 = Matrix::Zero(ncoords, ncoords);
  Matrix H_Pulay = Matrix::Zero(ncoords, ncoords);
  cout << "BEFORE one-body contributions to the hessian\n";
  //////////
  // one-body contributions to the hessian
  //////////
  auto T2 = compute_1body_ints_deriv<Operator::kinetic>(2, obs, atoms);
  auto V2 = compute_1body_ints_deriv<Operator::nuclear>(2, obs, atoms);
  for (auto row = 0, i = 0; row != ncoords; ++row) {
    for (auto col = row; col != ncoords; ++col, ++i) {
      auto hess = 2 * (T2[i] + V2[i]).cwiseProduct(D).sum();
      H1(row, col) += hess;
    }
  }

  cout << "BEFORE Pulay hessian\n";
  //////////
  // Pulay hessian
  //////////
  // orbital energy density
  //DiagonalMatrix evals_occalpha(evals.topRows(naocc)); // ORIGINAL CODE. COMMENTED.
  //DiagonalMatrix evals_occbeta(evals.topRows(nbocc)); // ORIGINAL CODE. COMMENTED.
  //Matrix W = Calpha_occ * evals_occalpha * Calpha_occ.transpose() + Cbeta_occ * evals_occbeta * Cbeta_occ.transpose();
  Matrix W = FockHWA + FockHWB;
  cout << "START printEigenmatrixCleanFmt(W)\n";
  printEigenmatrixCleanFmt(W);
  cout << "END printEigenmatrixCleanFmt(W)\n";
  auto S2 = compute_1body_ints_deriv<Operator::overlap>(2, obs, atoms);
  for (auto row = 0, i = 0; row != ncoords; ++row) {
    for (auto col = row; col != ncoords; ++col, ++i) {
      //auto hess = 2 * S2[i].cwiseProduct(W).sum(); // ORIGINAL CODE. COMMENTED. 2 SEEMS TO BE DUE TO SPIN.
      auto hess = S2[i].cwiseProduct(W).sum();
      H_Pulay(row, col) -= hess;
    }
  }

  cout << "** 1-body hessian = ";
  for (auto row = 0, i = 0; row != ncoords; ++row) {
    for (auto col = row; col != ncoords; ++col) {
      cout << H1(row, col) << " ";
    }
  }
  cout << std::endl;

  cout << "** Pulay hessian = ";
  for (auto row = 0, i = 0; row != ncoords; ++row) {
    for (auto col = row; col != ncoords; ++col) {
      cout << H_Pulay(row, col) << " ";
    }
  }
  cout << std::endl;
#endif  // LIBINT2_DERIV_ONEBODY_ORDER > 1

#if LIBINT2_DERIV_ERI_ORDER > 1
  cout << "BEFORE compute 2-e forces\n";
  // compute 2-e forces
  Matrix H2 = Matrix::Zero(ncoords, ncoords);

  //////////
  // two-body contributions to the forces
  //////////
  auto G2 = compute_2body_fock_deriv<2>(obs, atoms, D);
  for (auto row = 0, i = 0; row != ncoords; ++row) {
    for (auto col = row; col != ncoords; ++col, ++i) {
      // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
      auto hess = G2[i].cwiseProduct(D).sum();
      H2(row, col) += hess;
    }
  }

  cout << "** 2-body hessian = ";
  for (auto row = 0, i = 0; row != ncoords; ++row) {
    for (auto col = row; col != ncoords; ++col) {
      cout << H2(row, col) << " ";
    }
  }
  cout << std::endl;
#endif

// if support 1-e and 2-e 2nd derivatives compute nuclear repulsion hessian and
// the total hessian
#if LIBINT2_DERIV_ONEBODY_ORDER > 1 && LIBINT2_DERIV_ERI_ORDER > 1
  cout << "BEFORE compute nuclear repulsion hessian\n";
  // compute nuclear repulsion hessian
  // NB only the upper triangle is computed!!!
  Matrix HN = Matrix::Zero(ncoords, ncoords);

  //////////
  // nuclear repulsion contribution to the hessian
  //////////
  for (auto a1 = 1; a1 != atoms.size(); ++a1) {
    const auto& atom1 = atoms[a1];
    for (auto a2 = 0; a2 < a1; ++a2) {
      const auto& atom2 = atoms[a2];

      auto x12 = atom1.x - atom2.x;
      auto y12 = atom1.y - atom2.y;
      auto z12 = atom1.z - atom2.z;
      auto x12_2 = x12 * x12;
      auto y12_2 = y12 * y12;
      auto z12_2 = z12 * z12;
      auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
      auto r12 = sqrt(r12_2);
      auto r12_5 = r12 * r12_2 * r12_2;

      auto z1z2_over_r12_5 =
          atom1.atomic_number * atom2.atomic_number / r12_5;

      HN(3*a1 + 0, 3*a1 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
      HN(3*a1 + 1, 3*a1 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
      HN(3*a1 + 2, 3*a1 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
      HN(3*a1 + 0, 3*a1 + 1) += z1z2_over_r12_5 * (3*x12*y12);
      HN(3*a1 + 0, 3*a1 + 2) += z1z2_over_r12_5 * (3*x12*z12);
      HN(3*a1 + 1, 3*a1 + 2) += z1z2_over_r12_5 * (3*y12*z12);

      HN(3*a2 + 0, 3*a2 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
      HN(3*a2 + 1, 3*a2 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
      HN(3*a2 + 2, 3*a2 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
      HN(3*a2 + 0, 3*a2 + 1) += z1z2_over_r12_5 * (3*x12*y12);
      HN(3*a2 + 0, 3*a2 + 2) += z1z2_over_r12_5 * (3*x12*z12);
      HN(3*a2 + 1, 3*a2 + 2) += z1z2_over_r12_5 * (3*y12*z12);

      HN(3*a2 + 0, 3*a1 + 0) -= z1z2_over_r12_5 * (3*x12_2 - r12_2);
      HN(3*a2 + 1, 3*a1 + 1) -= z1z2_over_r12_5 * (3*y12_2 - r12_2);
      HN(3*a2 + 2, 3*a1 + 2) -= z1z2_over_r12_5 * (3*z12_2 - r12_2);
      HN(3*a2 + 1, 3*a1 + 0) -= z1z2_over_r12_5 * (3*y12*x12);
      HN(3*a2 + 2, 3*a1 + 0) -= z1z2_over_r12_5 * (3*z12*x12);
      HN(3*a2 + 2, 3*a1 + 1) -= z1z2_over_r12_5 * (3*z12*y12);
      HN(3*a2 + 0, 3*a1 + 1) -= z1z2_over_r12_5 * (3*x12*y12);
      HN(3*a2 + 0, 3*a1 + 2) -= z1z2_over_r12_5 * (3*x12*z12);
      HN(3*a2 + 1, 3*a1 + 2) -= z1z2_over_r12_5 * (3*y12*z12);
    }
  }

  cout << "** nuclear repulsion hessian = ";
  for (auto row = 0, i = 0; row != ncoords; ++row) {
    for (auto col = row; col != ncoords; ++col) {
      cout << HN(row, col) << " ";
    }
  }
  cout << std::endl;

  auto H = H1 + H_Pulay + H2 + HN;
  cout << "** Hartree-Fock hessian = ";
  for (auto row = 0, i = 0; row != ncoords; ++row) {
    for (auto col = row; col != ncoords; ++col) {
      cout << H(row, col) << " ";
    }
  }
  cout << std::endl;

  cout << "---------------------------------------------------------\n";
  cout << "-------------------- LAST ANTC CYCLE --------------------\n";
  cout << "------------------ END COMPUTE HESSIAN ------------------\n";
  cout << "---------------------------------------------------------\n";

#endif
}
*/
//-------------------------------------------------------------------------------------------
//-- END REMOVED FROM RUNANT BEFORE LAST CYCLE ANTC CALL, I PUTTED IT THERE, NOT ORIGINAL ---
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
//------------------------- END RUN ANT -----------------------------------------------------
//-------------------------------------------------------------------------------------------

int **HartreeFockClass::alloc_2d_int(int rows, int cols) {
    int *data = (int *)malloc(rows*cols*sizeof(int));
    int **array= (int **)malloc(rows*sizeof(int*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

double **HartreeFockClass::alloc_2d_double(int rows, int cols) {
	double *data = (double *)malloc(rows*cols*sizeof(double));
	double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

void HartreeFockClass::printDoubleMatrix(double * mat, int ROW, int COL)
//void HartreeFockClass::printDoubleMatrix(double mat[], int ROW, int COL)
{
    //int ROW = *mat->size();
    cout<<"\n Printing Matrix : \n";
    for(int i=0 ; i<=ROW-1 ; i++)
    {
        for(int j=0 ; j<=COL-1 ; j++)
            //cout<< *(*(mat+i)+j)<<" ";
        	cout<< mat[COL*i+j]<<" ";
        cout<<endl;
    }
    cout<<endl;
}

void HartreeFockClass::writeEigenMatrixToCSVfile(string name, Matrix matrix)
//void HartreeFockClass::writeEigenMatrixToCSVfile(const char * name, Matrix matrix)
{
	//const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
	const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, "  ", "\n");
    ofstream file(name.c_str());
    file << matrix.format(CSVFormat);
 }

void HartreeFockClass::printEigenmatrix(Matrix &m1){
	cout << "ENTER printEigenmatrix\n";
	std::string sep = "----------------------------------------\n";
	//Matrix3d m1;
	//m1 << 1.111111, 2, 3.33333, 4, 5, 6, 7, 8.888888, 9;
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
	Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
	Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
	cout << m1 << sep;
	cout << m1.format(CommaInitFmt) << "\n" << sep;
	cout << m1.format(CleanFmt) << "\n" << sep;
	cout << m1.format(OctaveFmt) << "\n" << sep;
	cout << m1.format(HeavyFmt) << "\n" << sep;
	cout << "EXIT printEigenmatrix\n";
}


void HartreeFockClass::printEigenmatrixCommaInitFmt(Matrix &m1){
	cout << "ENTER printEigenmatrixCommaInitFmt\n";
	std::string sep = "----------------------------------------\n";
	//Matrix3d m1;
	//m1 << 1.111111, 2, 3.33333, 4, 5, 6, 7, 8.888888, 9;
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
	//cout << m1 << sep;
	cout << sep << m1.format(CommaInitFmt) << "\n" << sep;
	cout << "EXIT printEigenmatrixCommaInitFmt\n";
}

void HartreeFockClass::printEigenmatrixCleanFmt(Matrix &m1){
	cout << "ENTER printEigenmatrixCleanFmt\n";
	std::string sep = "----------------------------------------\n";
	//Matrix3d m1;
	//m1 << 1.111111, 2, 3.33333, 4, 5, 6, 7, 8.888888, 9;
	Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	//cout << m1 << sep;
	cout << sep << m1.format(CleanFmt) << "\n" << sep;
	cout << "EXIT printEigenmatrixCleanFmt\n";
}

void HartreeFockClass::printEigenmatrixOctaveFmt(Matrix &m1){
	cout << "ENTER printEigenmatrixOctaveFmt\n";
	std::string sep = "----------------------------------------\n";
	//Matrix3d m1;
	//m1 << 1.111111, 2, 3.33333, 4, 5, 6, 7, 8.888888, 9;
	Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
	//cout << m1 << sep;
	cout << sep << m1.format(OctaveFmt) << "\n" << sep;
	cout << "EXIT printEigenmatrixOctaveFmt\n";
}

void HartreeFockClass::printEigenmatrixHeavyFmt(Matrix &m1){
	cout << "ENTER printEigenmatrixHeavyFmt\n";
	std::string sep = "----------------------------------------\n";
	//Matrix3d m1;
	//m1 << 1.111111, 2, 3.33333, 4, 5, 6, 7, 8.888888, 9;
	Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
	//cout << m1 << sep;
	cout << sep << m1.format(HeavyFmt) << "\n" << sep;
	cout << "EXIT printEigenmatrixHeavyFmt\n";
}

std::vector<Atom> HartreeFockClass::read_geometry(const std::string& filename) {
  cout << "Will read geometry from " << filename << std::endl;
  std::ifstream is(filename);
  if (not is.good()) {
    char errmsg[256] = "Could not open file ";
    strncpy(errmsg + 20, filename.c_str(), 235);
    errmsg[255] = '\0';
    throw std::runtime_error(errmsg);
  }

  // to prepare for MPI parallelization, we will read the entire file into a
  // string that can be
  // broadcast to everyone, then converted to an std::istringstream object that
  // can be used just like std::ifstream
  std::ostringstream oss;
  oss << is.rdbuf();
  // use ss.str() to get the entire contents of the file as an std::string
  // broadcast
  // then make an std::istringstream in each process
  std::istringstream iss(oss.str());

  // check the extension: if .xyz, assume the standard XYZ format, otherwise
  // throw an exception
  if (filename.rfind(".xyz") != std::string::npos)
    return libint2::read_dotxyz(iss);
  else
    throw "only .xyz files are accepted";
}

// computes Superposition-Of-Atomic-Densities guess for the molecular density
// matrix
// in minimal basis; occupies subshells by smearing electrons evenly over the
// orbitals
Matrix HartreeFockClass::compute_soad(const std::vector<Atom>& atoms) {
  // compute number of atomic orbitals
  size_t nao = 0;
  for (const auto& atom : atoms) {
    const auto Z = atom.atomic_number;
    nao += libint2::sto3g_num_ao(Z);
  }

  // compute the minimal basis density
  Matrix D = Matrix::Zero(nao, nao);
  size_t ao_offset = 0;  // first AO of this atom
  for (const auto& atom : atoms) {
    const auto Z = atom.atomic_number;
    const auto& occvec = libint2::sto3g_ao_occupation_vector(Z);
    for(const auto& occ: occvec) {
      D(ao_offset, ao_offset) = occ;
      ++ao_offset;
    }
  }

  return D * 0.5;  // we use densities normalized to # of electrons/2
}

Matrix HartreeFockClass::compute_shellblock_norm(const BasisSet& obs, const Matrix& A) {
  const auto nsh = obs.size();
  Matrix Ash(nsh, nsh);

  auto shell2bf = obs.shell2bf();
  for (size_t s1 = 0; s1 != nsh; ++s1) {
    const auto& s1_first = shell2bf[s1];
    const auto& s1_size = obs[s1].size();
    for (size_t s2 = 0; s2 != nsh; ++s2) {
      const auto& s2_first = shell2bf[s2];
      const auto& s2_size = obs[s2].size();

      Ash(s1, s2) = A.block(s1_first, s2_first, s1_size, s2_size)
                        .lpNorm<Eigen::Infinity>();
    }
  }

  return Ash;
}

template <Operator obtype>
std::array<Matrix, libint2::operator_traits<obtype>::nopers> HartreeFockClass::compute_1body_ints(
    const BasisSet& obs, const std::vector<Atom>& atoms) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  using libint2::nthreads;
  typedef std::array<Matrix, libint2::operator_traits<obtype>::nopers>
      result_type;
  const unsigned int nopers = libint2::operator_traits<obtype>::nopers;
  result_type result;
  for (auto& r : result) r = Matrix::Zero(n, n);

  // construct the 1-body integrals engine
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(obtype, obs.max_nprim(), obs.max_l(), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical
  // charges
  if (obtype == Operator::nuclear) {
    engines[0].set_params(libint2::make_point_charges(atoms));
  }
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = obs.shell2bf();

  auto compute = [&](int thread_id) {

    const auto& buf = engines[thread_id].results();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over
    // Hermitian operators: (1|2) = (2|1)
    for (auto s1 = 0l, s12 = 0l; s1 != nshells; ++s1) {
      auto bf1 = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();

      auto s1_offset = s1 * (s1+1) / 2;
      for (auto s2: obs_shellpair_list[s1]) {
        auto s12 = s1_offset + s2;
        if (s12 % nthreads != thread_id) continue;

        auto bf2 = shell2bf[s2];
        auto n2 = obs[s2].size();

        auto n12 = n1 * n2;

        // compute shell pair; return is the pointer to the buffer
        engines[thread_id].compute(obs[s1], obs[s2]);

        for (unsigned int op = 0; op != nopers; ++op) {
          // "map" buffer to a const Eigen Matrix, and copy it to the
          // corresponding blocks of the result
          Eigen::Map<const Matrix> buf_mat(buf[op], n1, n2);
          result[op].block(bf1, bf2, n1, n2) = buf_mat;
          if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding
                         // {s2,s1} block, note the transpose!
            result[op].block(bf2, bf1, n2, n1) = buf_mat.transpose();
        }
      }
    }
  };  // compute lambda

  libint2::parallel_do(compute);

  return result;
}

#if LIBINT2_DERIV_ONEBODY_ORDER
template <Operator obtype>
std::vector<Matrix> HartreeFockClass::compute_1body_ints_deriv(unsigned deriv_order,
                                             const BasisSet& obs,
                                             const std::vector<Atom>& atoms) {
  using libint2::nthreads;
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  constexpr auto nopers = libint2::operator_traits<obtype>::nopers;
  const auto nresults =
      nopers * libint2::num_geometrical_derivatives(atoms.size(), deriv_order);
  typedef std::vector<Matrix> result_type;
  result_type result(nresults);
  for (auto& r : result) r = Matrix::Zero(n, n);

  // construct the 1-body integrals engine
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] =
      libint2::Engine(obtype, obs.max_nprim(), obs.max_l(), deriv_order);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical
  // charges
  if (obtype == Operator::nuclear) {
    std::vector<std::pair<double, std::array<double, 3>>> q;
    for (const auto& atom : atoms) {
      q.push_back({static_cast<double>(atom.atomic_number),
                   {{atom.x, atom.y, atom.z}}});
    }
    engines[0].set_params(q);
  }
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = obs.shell2bf();
  auto shell2atom = obs.shell2atom(atoms);

  const auto natoms = atoms.size();
  const auto two_times_ncoords = 6*natoms;
  const auto nderivcenters_shset =
      2 + ((obtype == Operator::nuclear) ? natoms : 0);

  auto compute = [&](int thread_id) {

    const auto& buf = engines[thread_id].results();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over
    // Hermitian operators: (1|2) = (2|1)
    for (auto s1 = 0l, s12 = 0l; s1 != nshells; ++s1) {
      auto bf1 = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();
      auto atom1 = shell2atom[s1];
      assert(atom1 != -1);

      auto s1_offset = s1 * (s1+1) / 2;
      for (auto s2: obs_shellpair_list[s1]) {
        auto s12 = s1_offset + s2;
        if (s12 % nthreads != thread_id) continue;

        auto bf2 = shell2bf[s2];
        auto n2 = obs[s2].size();
        auto atom2 = shell2atom[s2];

        auto n12 = n1 * n2;

        // compute shell pair; return is the pointer to the buffer
        engines[thread_id].compute(obs[s1], obs[s2]);

        // "copy" lambda copies shell set \c idx to the operator matrix with
        // index \c op
        auto add_shellset_to_dest = [&](std::size_t op, std::size_t idx,
                                        double scale = 1.0) {
          // "map" buffer to a const Eigen Matrix, and copy it to the
          // corresponding blocks of the result
          Eigen::Map<const Matrix> buf_mat(buf[idx], n1, n2);
          if (scale == 1.0)
            result[op].block(bf1, bf2, n1, n2) += buf_mat;
          else
            result[op].block(bf1, bf2, n1, n2) += scale * buf_mat;
          if (s1 != s2) {  // if s1 >= s2, copy {s1,s2} to the corresponding
                           // {s2,s1} block, note the transpose!
            if (scale == 1.0)
              result[op].block(bf2, bf1, n2, n1) += buf_mat.transpose();
            else
              result[op].block(bf2, bf1, n2, n1) += scale * buf_mat.transpose();
          }
        };

        switch (deriv_order) {
          case 0:
            for (std::size_t op = 0; op != nopers; ++op) {
              add_shellset_to_dest(op, op);
            }
            break;

          // map deriv quanta for this shell pair to the overall deriv quanta
          //
          // easiest to explain with example:
          // in sto-3g water shells 0 1 2 sit on atom 0, shells 3 and 4 on atoms
          // 1 and 2 respectively
          // each call to engine::compute for nuclear ints will return
          // derivatives
          // with respect to 15 coordinates, obtained as 3 (x,y,z) times 2 + 3 =
          // 5 centers
          // (2 centers on which shells sit + 3 nuclear charges)
          // (for overlap, kinetic, and emultipole ints we there are only 6
          // coordinates
          //  since the operator is coordinate-independent, or derivatives with
          //  respect to
          //  the operator coordinates are not computed)
          //

          case 1: {
            std::size_t shellset_idx = 0;
            for (auto c = 0; c != nderivcenters_shset; ++c) {
              auto atom = (c == 0) ? atom1 : ((c == 1) ? atom2 : c - 2);
              auto op_start = 3 * atom * nopers;
              auto op_fence = op_start + nopers;
              for (auto xyz = 0; xyz != 3;
                   ++xyz, op_start += nopers, op_fence += nopers) {
                for (unsigned int op = op_start; op != op_fence;
                     ++op, ++shellset_idx) {
                  add_shellset_to_dest(op, shellset_idx);
                }
              }
            }
          } break;

          case 2: {
            //
            // must pay attention to symmetry when computing 2nd and higher-order derivs
            // e.g. d2 (s1|s2) / dX dY involves several cases:
            // 1. only s1 (or only s2) depends on X AND Y (i.e. X and Y refer to same atom) =>
            //    d2 (s1|s2) / dX dY = (d2 s1 / dX dY | s2)
            // 2. s1 depends on X only, s2 depends on Y only (or vice versa) =>
            //    d2 (s1|s2) / dX dY = (d s1 / dX | d s2 / dY)
            // 3. s1 AND s2 depend on X AND Y (i.e. X and Y refer to same atom) =>
            //    case A: X != Y
            //    d2 (s1|s2) / dX dY = (d2 s1 / dX dY | s2) + (d s1 / dX | d s2 / dY)
            //      + (d s1 / dY | d s2 / dX) + (s1| d2 s2 / dX dY )
            //    case B: X == Y
            //    d2 (s1|s2) / dX2 = (d2 s1 / dX2 | s2) + 2 (d s1 / dX | d s2 / dX)
            //      + (s1| d2 s2 / dX2 )

            // computes upper triangle index
            // n2 = matrix size times 2
            // i,j = (unordered) indices
#define upper_triangle_index(n2, i, j) \
  (std::min((i), (j))) * ((n2) - (std::min((i), (j))) - 1) / 2 + \
            (std::max((i), (j)))

            // look over shellsets in the order in which they appear
            std::size_t shellset_idx = 0;
            for (auto c1 = 0; c1 != nderivcenters_shset; ++c1) {
              auto a1 = (c1 == 0) ? atom1 : ((c1 == 1) ? atom2 : c1 - 2);
              auto coord1 = 3 * a1;
              for (auto xyz1 = 0; xyz1 != 3; ++xyz1, ++coord1) {

                for (auto c2 = c1; c2 != nderivcenters_shset; ++c2) {
                  auto a2 = (c2 == 0) ? atom1 : ((c2 == 1) ? atom2 : c2 - 2);
                  auto xyz2_start = (c1 == c2) ? xyz1 : 0;
                  auto coord2 = 3 * a2 + xyz2_start;
                  for (auto xyz2 = xyz2_start; xyz2 != 3;
                       ++xyz2, ++coord2) {

                    double scale = (coord1 == coord2 && c1 != c2) ? 2.0 : 1.0;

                    const auto coord12 =
                        upper_triangle_index(two_times_ncoords, coord1, coord2);
                    auto op_start = coord12 * nopers;
                    auto op_fence = op_start + nopers;
                    for (auto op = op_start; op != op_fence;
                        ++op, ++shellset_idx) {
                      add_shellset_to_dest(op, shellset_idx, scale);
                    }
                  }
                }
              }
            }
          } break;
#undef upper_triangle_index

          default: {
            assert(false && "not yet implemented");

            using ShellSetDerivIterator =
                libint2::FixedOrderedIntegerPartitionIterator<
                    std::vector<unsigned int>>;
            ShellSetDerivIterator shellset_diter(deriv_order,
                                                 nderivcenters_shset);
            while (shellset_diter) {
              const auto& deriv = *shellset_diter;
            }
          }
        }  // copy shell block switch

      }  // s2 <= s1
    }    // s1
  };  // compute lambda

  libint2::parallel_do(compute);

  return result;
}
#endif

template <libint2::Operator Kernel>
Matrix HartreeFockClass::compute_schwarz_ints(
    const BasisSet& bs1, const BasisSet& _bs2, bool use_2norm,
    typename libint2::operator_traits<Kernel>::oper_params_type params) {
	cout << "ENTER compute_schwarz_ints\n";
  const BasisSet& bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  Matrix K = Matrix::Zero(nsh1, nsh2);

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  using libint2::nthreads;
  std::vector<Engine> engines(nthreads);

  // !!! very important: cannot screen primitives in Schwarz computation !!!
  auto epsilon = 0.;
  engines[0] = Engine(Kernel, std::max(bs1.max_nprim(), bs2.max_nprim()),
                      std::max(bs1.max_l(), bs2.max_l()), 0, epsilon, params);
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  cout << "computing Schwarz bound prerequisites (kernel=" << (int)Kernel
            << ") ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  auto compute = [&](int thread_id) {

    const auto& buf = engines[thread_id].results();

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      auto n1 = bs1[s1].size();  // number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
        if (s12 % nthreads != thread_id) continue;

        auto n2 = bs2[s2].size();
        auto n12 = n1 * n2;

        engines[thread_id].compute2<Kernel, BraKet::xx_xx, 0>(bs1[s1], bs2[s2],
                                                              bs1[s1], bs2[s2]);
        assert(buf[0] != nullptr &&
               "to compute Schwarz ints turn off primitive screening");

        // to apply Schwarz inequality to individual integrals must use the diagonal elements
        // to apply it to sets of functions (e.g. shells) use the whole shell-set of ints here
        Eigen::Map<const Matrix> buf_mat(buf[0], n12, n12);
        auto norm2 = use_2norm ? buf_mat.norm()
                               : buf_mat.lpNorm<Eigen::Infinity>();
        K(s1, s2) = std::sqrt(norm2);
        if (bs1_equiv_bs2) K(s2, s1) = K(s1, s2);
      }
    }
  };  // thread lambda

  libint2::parallel_do(compute);

  timer.stop(0);
  cout << "done (" << timer.read(0) << " s)" << std::endl;

  cout << "EXIT compute_schwarz_ints\n";

  return K;
}

Matrix HartreeFockClass::compute_do_ints(const BasisSet& bs1, const BasisSet& bs2,
                       bool use_2norm) {
  return compute_schwarz_ints<libint2::Operator::delta>(bs1, bs2, use_2norm);
}

shellpair_list_t HartreeFockClass::compute_shellpair_list(const BasisSet& bs1,
                                        const BasisSet& _bs2,
                                        const double threshold) {
  const BasisSet& bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  using libint2::nthreads;

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  std::vector<Engine> engines;
  engines.reserve(nthreads);
  engines.emplace_back(Operator::overlap,
                       std::max(bs1.max_nprim(), bs2.max_nprim()),
                       std::max(bs1.max_l(), bs2.max_l()), 0);
  for (size_t i = 1; i != nthreads; ++i) {
    engines.push_back(engines[0]);
  }

  cout << "computing non-negligible shell-pair list ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  shellpair_list_t result;

  std::mutex mx;

  auto compute = [&](int thread_id) {

    auto& engine = engines[thread_id];
    const auto& buf = engine.results();

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      mx.lock();
      if (result.find(s1) == result.end())
        result.insert(std::make_pair(s1, std::vector<size_t>()));
      mx.unlock();

      auto n1 = bs1[s1].size();  // number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
        if (s12 % nthreads != thread_id) continue;

        auto on_same_center = (bs1[s1].O == bs2[s2].O);
        bool significant = on_same_center;
        if (not on_same_center) {
          auto n2 = bs2[s2].size();
          engines[thread_id].compute(bs1[s1], bs2[s2]);
          Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
          auto norm = buf_mat.norm();
          significant = (norm >= threshold);
        }

        if (significant) {
          mx.lock();
          result[s1].emplace_back(s2);
          mx.unlock();
        }
      }
    }
  };  // end of compute

  libint2::parallel_do(compute);

  // resort shell list in increasing order, i.e. result[s][s1] < result[s][s2]
  // if s1 < s2
  auto sort = [&](int thread_id) {
    for (auto s1 = 0l; s1 != nsh1; ++s1) {
      if (s1 % nthreads == thread_id) {
        auto& list = result[s1];
        std::sort(list.begin(), list.end());
      }
    }
  };  // end of sort

  libint2::parallel_do(sort);

  timer.stop(0);
  cout << "done (" << timer.read(0) << " s)" << std::endl;

  return result;
}

// returns {X,X^{-1},rank,A_condition_number,result_A_condition_number}, where
// X is the generalized square-root-inverse such that X.transpose() * A * X = I
//
// if symmetric is true, produce "symmetric" sqrtinv: X = U.A_evals_sqrtinv .
// U.transpose()),
// else produce "canonical" sqrtinv: X = U.A_evals_sqrtinv
// where U are eigenvectors of A
// rows and cols of symmetric X are equivalent; for canonical X the rows are
// original basis (AO),
// cols are transformed basis ("orthogonal" AO)
//
// A is conditioned to max_condition_number

//std::tuple<Matrix, Matrix, size_t, double, double> HartreeFockClass::gensqrtinv(
//    const Matrix& S, bool symmetric = false,
//    double max_condition_number = 1e8) {
std::tuple<Matrix, Matrix, size_t, double, double> HartreeFockClass::gensqrtinv(
    const Matrix& S, bool symmetric,
    double max_condition_number ) {
  Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(S);
  auto U = eig_solver.eigenvectors();
  auto s = eig_solver.eigenvalues();
  auto s_max = s.maxCoeff();
  auto condition_number = std::min(
      s_max / std::max(s.minCoeff(), std::numeric_limits<double>::min()),
      1.0 / std::numeric_limits<double>::epsilon());
  auto threshold = s_max / max_condition_number;
  long n = s.rows();
  long n_cond = 0;
  for (long i = n - 1; i >= 0; --i) {
    if (s(i) >= threshold) {
      ++n_cond;
    } else
      i = 0;  // skip rest since eigenvalues are in ascending order
  }

  auto sigma = s.bottomRows(n_cond);
  auto result_condition_number = sigma.maxCoeff() / sigma.minCoeff();
  auto sigma_sqrt = sigma.array().sqrt().matrix().asDiagonal();
  auto sigma_invsqrt = sigma.array().sqrt().inverse().matrix().asDiagonal();

  // make canonical X/Xinv
  auto U_cond = U.block(0, n - n_cond, n, n_cond);
  Matrix X = U_cond * sigma_invsqrt;
  Matrix Xinv = U_cond * sigma_sqrt;
  // convert to symmetric, if needed
  if (symmetric) {
    X = X * U_cond.transpose();
    Xinv = Xinv * U_cond.transpose();
  }
  return std::make_tuple(X, Xinv, size_t(n_cond), condition_number,
                         result_condition_number);
}

std::tuple<Matrix, Matrix, double> HartreeFockClass::conditioning_orthogonalizer(
    const Matrix& S, double S_condition_number_threshold) {
  size_t obs_rank;
  double S_condition_number;
  double XtX_condition_number;
  Matrix X, Xinv;

  assert(S.rows() == S.cols());

  std::tie(X, Xinv, obs_rank, S_condition_number, XtX_condition_number) =
      gensqrtinv(S, false, S_condition_number_threshold);
  auto obs_nbf_omitted = (long)S.rows() - (long)obs_rank;
  cout << "overlap condition number = " << S_condition_number;
  if (obs_nbf_omitted > 0)
    cout << " (dropped " << obs_nbf_omitted << " "
              << (obs_nbf_omitted > 1 ? "fns" : "fn") << " to reduce to "
              << XtX_condition_number << ")";
  cout << std::endl;

  if (obs_nbf_omitted > 0) {
    Matrix should_be_I = X.transpose() * S * X;
    Matrix I = Matrix::Identity(should_be_I.rows(), should_be_I.cols());
    cout << "||X^t * S * X - I||_2 = " << (should_be_I - I).norm()
              << " (should be 0)" << std::endl;
  }

  return std::make_tuple(X, Xinv, XtX_condition_number);
}

Matrix HartreeFockClass::compute_2body_2index_ints(const BasisSet& bs) {
  using libint2::nthreads;
  const auto n = bs.nbf();
  const auto nshells = bs.size();
  Matrix result = Matrix::Zero(n, n);

  // build engines for each thread
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] =
      Engine(libint2::Operator::coulomb, bs.max_nprim(), bs.max_l(), 0);
  engines[0].set_braket(BraKet::xs_xs);
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = bs.shell2bf();
  auto unitshell = Shell::unit();

  auto compute = [&](int thread_id) {

    auto& engine = engines[thread_id];
    const auto& buf = engine.results();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over
    // Hermitian operators: (1|2) = (2|1)
    for (auto s1 = 0l, s12 = 0l; s1 != nshells; ++s1) {
      auto bf1 = shell2bf[s1];  // first basis function in this shell
      auto n1 = bs[s1].size();

      for (auto s2 = 0; s2 <= s1; ++s2, ++s12) {
        if (s12 % nthreads != thread_id) continue;

        auto bf2 = shell2bf[s2];
        auto n2 = bs[s2].size();

        // compute shell pair; return is the pointer to the buffer
        engine.compute(bs[s1], bs[s2]);
        if (buf[0] == nullptr)
          continue; // if all integrals screened out, skip to next shell set

        // "map" buffer to a const Eigen Matrix, and copy it to the
        // corresponding blocks of the result
        Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
        result.block(bf1, bf2, n1, n2) = buf_mat;
        if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                       // block, note the transpose!
          result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
      }
    }
  };  // compute lambda

  libint2::parallel_do(compute);

  return result;
}

Matrix HartreeFockClass::compute_2body_fock(const BasisSet& obs, const Matrix& D,
                          double precision, const Matrix& Schwarz) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  using libint2::nthreads;
  std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));

  const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
  Matrix D_shblk_norm =
      compute_shellblock_norm(obs, D);  // matrix of infty-norms of shell blocks

  auto fock_precision = precision;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  auto max_nprim = obs.max_nprim();
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
  auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  cout << "compute_2body_fock:precision = " << precision << std::endl;
  cout << "Engine::precision = " << engines[0].precision() << std::endl;
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = obs.shell2bf();

  auto lambda = [&](int thread_id) {

    auto& engine = engines[thread_id];
    auto& g = G[thread_id];
    const auto& buf = engine.results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto& timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell

      for (const auto& s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();

          const auto Dnorm123 =
              do_schwarz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto& s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break;  // for each s3, s4 are stored in monotonically increasing
                      // order

            if ((s1234++) % nthreads != thread_id) continue;

            const auto Dnorm1234 =
                do_schwarz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwarz_screen &&
                Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                obs[s1], obs[s2], obs[s3], obs[s4]);
            const auto* buf_1234 = buf[0];
            if (buf_1234 == nullptr)
              continue; // if all integrals screened out, skip to next quartet

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            // 1) each shell set of integrals contributes up to 6 shell sets of
            // the Fock matrix:
            // LOOKS LIKE THE NEGATIVE TERMS ARE EXCHANGE TERMS AND...
            // THE POSITIVE TERMS ARE COULOMB TERMS.
            //    F(a,b) += (ab|cd) * D(c,d)
            //    F(c,d) += (ab|cd) * D(a,b)
            //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
            //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
            //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
            //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
            // 2) each permutationally-unique integral (shell set) must be
            // scaled by its degeneracy,
            //    i.e. the number of the integrals/sets equivalent to it
            // 3) the end result must be symmetrized
            for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for (auto f3 = 0; f3 != n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf_1234[f1234];

                    const auto value_scal_by_deg = value * s1234_deg;

                    g(bf1, bf2) += D(bf3, bf4) * value_scal_by_deg;
                    g(bf3, bf4) += D(bf1, bf2) * value_scal_by_deg;
                    g(bf1, bf3) -= 0.25 * D(bf2, bf4) * value_scal_by_deg;
                    g(bf2, bf4) -= 0.25 * D(bf1, bf3) * value_scal_by_deg;
                    g(bf1, bf4) -= 0.25 * D(bf2, bf3) * value_scal_by_deg;
                    g(bf2, bf3) -= 0.25 * D(bf1, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }

  };  // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    G[0] += G[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto& t : timers) {
    time_for_ints += t.read(0);
  }
  cout << "time for integrals = " << time_for_ints << std::endl;
  for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
#endif

  Matrix GG = 0.5 * (G[0] + G[0].transpose());

  cout << "# of integrals = " << num_ints_computed << std::endl;

  // symmetrize the result and return
  return GG;
}

VectorOfMatrices HartreeFockClass::compute_2body_fock_uhf(const BasisSet& obs, const Matrix& Dalpha, const Matrix& Dbeta,
                          double precision, const Matrix& Schwarz) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  using libint2::nthreads;
  std::vector<Matrix> Galpha(nthreads, Matrix::Zero(n, n));
  std::vector<Matrix> Gbeta(nthreads, Matrix::Zero(n, n));

  const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
  Matrix Dalpha_shblk_norm =
      compute_shellblock_norm(obs, Dalpha);  // matrix of infty-norms of shell blocks
  Matrix Dbeta_shblk_norm =
      compute_shellblock_norm(obs, Dbeta);  // matrix of infty-norms of shell blocks

  auto fock_precision = precision;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  auto max_nprim = obs.max_nprim();
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
  auto engine_precision = std::min(fock_precision / (std::max(Dalpha_shblk_norm.maxCoeff(),Dbeta_shblk_norm.maxCoeff())),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  cout << "compute_2body_fock:precision = " << precision << std::endl;
  cout << "Engine::precision = " << engines[0].precision() << std::endl;
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = obs.shell2bf();

  auto lambda = [&](int thread_id) {

    auto& engine = engines[thread_id];
    auto& galpha = Galpha[thread_id];
    auto& gbeta = Gbeta[thread_id];
    const auto& buf = engine.results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto& timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell

      for (const auto& s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        const auto Dalphanorm12 = do_schwarz_screen ? Dalpha_shblk_norm(s1, s2) : 0.;
        const auto Dbetanorm12 = do_schwarz_screen ? Dbeta_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();

          const auto Dalphanorm123 =
              do_schwarz_screen
                  ? std::max(Dalpha_shblk_norm(s1, s3),
                             std::max(Dalpha_shblk_norm(s2, s3), Dalphanorm12))
                  : 0.;
          const auto Dbetanorm123 =
              do_schwarz_screen
                  ? std::max(Dbeta_shblk_norm(s1, s3),
                             std::max(Dbeta_shblk_norm(s2, s3), Dbetanorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto& s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break;  // for each s3, s4 are stored in monotonically increasing
                      // order

            if ((s1234++) % nthreads != thread_id) continue;

            const auto Dalphanorm1234 =
                do_schwarz_screen
                    ? std::max(
                          Dalpha_shblk_norm(s1, s4),
                          std::max(Dalpha_shblk_norm(s2, s4),
                                   std::max(Dalpha_shblk_norm(s3, s4), Dalphanorm123)))
                    : 0.;
            const auto Dbetanorm1234 =
                do_schwarz_screen
                    ? std::max(
                          Dbeta_shblk_norm(s1, s4),
                          std::max(Dbeta_shblk_norm(s2, s4),
                                   std::max(Dbeta_shblk_norm(s3, s4), Dbetanorm123)))
                    : 0.;

            if (do_schwarz_screen &&
                Dalphanorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                obs[s1], obs[s2], obs[s3], obs[s4]);
            const auto* buf_1234 = buf[0];
            if (buf_1234 == nullptr)
              continue; // if all integrals screened out, skip to next quartet

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            // 1) each shell set of integrals contributes up to 6 shell sets of
            // the Fock matrix:
            //    F(a,b) += (ab|cd) * D(c,d)
            //    F(c,d) += (ab|cd) * D(a,b)
            //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
            //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
            //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
            //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
            // 2) each permutationally-unique integral (shell set) must be
            // scaled by its degeneracy,
            //    i.e. the number of the integrals/sets equivalent to it
            // 3) the end result must be symmetrized
            for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for (auto f3 = 0; f3 != n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf_1234[f1234];

                    const auto value_scal_by_deg = value * s1234_deg;

                    galpha(bf1, bf2) += (Dalpha(bf3, bf4) + Dbeta(bf3, bf4)) * value_scal_by_deg;
                    galpha(bf3, bf4) += (Dalpha(bf1, bf2) + Dbeta(bf1, bf2)) * value_scal_by_deg;
                    galpha(bf1, bf3) -= 0.25 * Dalpha(bf2, bf4) * value_scal_by_deg;
                    galpha(bf2, bf4) -= 0.25 * Dalpha(bf1, bf3) * value_scal_by_deg;
                    galpha(bf1, bf4) -= 0.25 * Dalpha(bf2, bf3) * value_scal_by_deg;
                    galpha(bf2, bf3) -= 0.25 * Dalpha(bf1, bf4) * value_scal_by_deg;

                    gbeta(bf1, bf2) += (Dalpha(bf3, bf4) + Dbeta(bf3, bf4)) * value_scal_by_deg;
                    gbeta(bf3, bf4) += (Dalpha(bf1, bf2) + Dbeta(bf1, bf2)) * value_scal_by_deg;
                    gbeta(bf1, bf3) -= 0.25 * Dbeta(bf2, bf4) * value_scal_by_deg;
                    gbeta(bf2, bf4) -= 0.25 * Dbeta(bf1, bf3) * value_scal_by_deg;
                    gbeta(bf1, bf4) -= 0.25 * Dbeta(bf2, bf3) * value_scal_by_deg;
                    gbeta(bf2, bf3) -= 0.25 * Dbeta(bf1, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }

  };  // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    Galpha[0] += Galpha[i];
    Gbeta[0] += Gbeta[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto& t : timers) {
    time_for_ints += t.read(0);
  }
  cout << "time for integrals = " << time_for_ints << std::endl;
  for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
#endif

  // SYMMETRIZE RESULT.
  Matrix GGalpha = 0.5 * (Galpha[0] + Galpha[0].transpose());
  Matrix GGbeta = 0.5 * (Gbeta[0] + Gbeta[0].transpose());

  cout << "# of integrals = " << num_ints_computed << std::endl;

  VectorOfMatrices GGvectorAB(2);
  GGvectorAB = {GGalpha,GGbeta};
  // symmetrize the result and return
  return GGvectorAB;
}



VectorOfMatrices HartreeFockClass::compute_2body_JK_uhf(const BasisSet& obs, const Matrix& Dalpha, const Matrix& Dbeta,
                          double precision, const Matrix& Schwarz) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  using libint2::nthreads;
  std::vector<Matrix> J(nthreads, Matrix::Zero(n, n));
  std::vector<Matrix> Kalpha(nthreads, Matrix::Zero(n, n));
  std::vector<Matrix> Kbeta(nthreads, Matrix::Zero(n, n));

  const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
  Matrix Dalpha_shblk_norm =
      compute_shellblock_norm(obs, Dalpha);  // matrix of infty-norms of shell blocks
  Matrix Dbeta_shblk_norm =
      compute_shellblock_norm(obs, Dbeta);  // matrix of infty-norms of shell blocks

  auto fock_precision = precision;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  auto max_nprim = obs.max_nprim();
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
  auto engine_precision = std::min(fock_precision / (std::max(Dalpha_shblk_norm.maxCoeff(),Dbeta_shblk_norm.maxCoeff())),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  cout << "compute_2body_fock:precision = " << precision << std::endl;
  cout << "Engine::precision = " << engines[0].precision() << std::endl;
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = obs.shell2bf();

  auto lambda = [&](int thread_id) {

    auto& engine = engines[thread_id];
    auto& j = J[thread_id];
    auto& kalpha = Kalpha[thread_id];
    auto& kbeta = Kbeta[thread_id];
    const auto& buf = engine.results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto& timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell

      for (const auto& s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        const auto Dalphanorm12 = do_schwarz_screen ? Dalpha_shblk_norm(s1, s2) : 0.;
        const auto Dbetanorm12 = do_schwarz_screen ? Dbeta_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();

          const auto Dalphanorm123 =
              do_schwarz_screen
                  ? std::max(Dalpha_shblk_norm(s1, s3),
                             std::max(Dalpha_shblk_norm(s2, s3), Dalphanorm12))
                  : 0.;
          const auto Dbetanorm123 =
              do_schwarz_screen
                  ? std::max(Dbeta_shblk_norm(s1, s3),
                             std::max(Dbeta_shblk_norm(s2, s3), Dbetanorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto& s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break;  // for each s3, s4 are stored in monotonically increasing
                      // order

            if ((s1234++) % nthreads != thread_id) continue;

            const auto Dalphanorm1234 =
                do_schwarz_screen
                    ? std::max(
                          Dalpha_shblk_norm(s1, s4),
                          std::max(Dalpha_shblk_norm(s2, s4),
                                   std::max(Dalpha_shblk_norm(s3, s4), Dalphanorm123)))
                    : 0.;
            const auto Dbetanorm1234 =
                do_schwarz_screen
                    ? std::max(
                          Dbeta_shblk_norm(s1, s4),
                          std::max(Dbeta_shblk_norm(s2, s4),
                                   std::max(Dbeta_shblk_norm(s3, s4), Dbetanorm123)))
                    : 0.;

            if (do_schwarz_screen &&
                Dalphanorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                obs[s1], obs[s2], obs[s3], obs[s4]);
            const auto* buf_1234 = buf[0];
            if (buf_1234 == nullptr)
              continue; // if all integrals screened out, skip to next quartet

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            // 1) each shell set of integrals contributes up to 6 shell sets of
            // the Fock matrix:
            //    F(a,b) += (ab|cd) * D(c,d)
            //    F(c,d) += (ab|cd) * D(a,b)
            //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
            //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
            //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
            //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
            // 2) each permutationally-unique integral (shell set) must be
            // scaled by its degeneracy,
            //    i.e. the number of the integrals/sets equivalent to it
            // 3) the end result must be symmetrized
            for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for (auto f3 = 0; f3 != n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf_1234[f1234];

                    const auto value_scal_by_deg = value * s1234_deg;

                    j(bf1, bf2) += (Dalpha(bf3, bf4) + Dbeta(bf3, bf4)) * value_scal_by_deg;
                    j(bf3, bf4) += (Dalpha(bf1, bf2) + Dbeta(bf1, bf2)) * value_scal_by_deg;
                    //Jalpha(bf1, bf2) += (Dalpha(bf3, bf4) + Dbeta(bf3, bf4)) * value_scal_by_deg;
                    //Jalpha(bf3, bf4) += (Dalpha(bf1, bf2) + Dbeta(bf1, bf2)) * value_scal_by_deg;
                    kalpha(bf1, bf3) -= 0.25 * Dalpha(bf2, bf4) * value_scal_by_deg;
                    kalpha(bf2, bf4) -= 0.25 * Dalpha(bf1, bf3) * value_scal_by_deg;
                    kalpha(bf1, bf4) -= 0.25 * Dalpha(bf2, bf3) * value_scal_by_deg;
                    kalpha(bf2, bf3) -= 0.25 * Dalpha(bf1, bf4) * value_scal_by_deg;

                    //Jbeta(bf1, bf2) += (Dalpha(bf3, bf4) + Dbeta(bf3, bf4)) * value_scal_by_deg;
                    //Jbeta(bf3, bf4) += (Dalpha(bf1, bf2) + Dbeta(bf1, bf2)) * value_scal_by_deg;
                    kbeta(bf1, bf3) -= 0.25 * Dbeta(bf2, bf4) * value_scal_by_deg;
                    kbeta(bf2, bf4) -= 0.25 * Dbeta(bf1, bf3) * value_scal_by_deg;
                    kbeta(bf1, bf4) -= 0.25 * Dbeta(bf2, bf3) * value_scal_by_deg;
                    kbeta(bf2, bf3) -= 0.25 * Dbeta(bf1, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }

  };  // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
	J[0] += J[i];
    Kalpha[0] += Kalpha[i];
    Kbeta[0] += Kbeta[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto& t : timers) {
    time_for_ints += t.read(0);
  }
  cout << "time for integrals = " << time_for_ints << std::endl;
  for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
#endif

  // SYMMETRIZE RESULT.
  Matrix JJ = 0.5 * (J[0] + J[0].transpose());
  Matrix KKalpha = 0.5 * (Kalpha[0] + Kalpha[0].transpose());
  Matrix KKbeta = 0.5 * (Kbeta[0] + Kbeta[0].transpose());

  cout << "# of integrals = " << num_ints_computed << std::endl;

  VectorOfMatrices JKvectorAB(3);
  JKvectorAB = {JJ,KKalpha,KKbeta};
  // symmetrize the result and return
  return JKvectorAB;
}




#if LIBINT2_DERIV_ERI_ORDER
template <unsigned deriv_order>
std::vector<Matrix> HartreeFockClass::compute_2body_fock_deriv(const BasisSet& obs,
                                             const std::vector<Atom>& atoms,
                                             const Matrix& D, double precision,
                                             const Matrix& Schwarz) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  const auto nderiv_shellset =
      libint2::num_geometrical_derivatives(4, deriv_order); // # of derivs for each shell quartet
  const auto nderiv = libint2::num_geometrical_derivatives(
      atoms.size(), deriv_order);  // total # of derivs
  const auto ncoords_times_two = (atoms.size() * 3) * 2;
  using libint2::nthreads;
  std::vector<Matrix> G(nthreads * nderiv, Matrix::Zero(n, n));

  const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
  Matrix D_shblk_norm =
      compute_shellblock_norm(obs, D);  // matrix of infty-norms of shell blocks

  auto fock_precision = precision;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  auto max_nprim = obs.max_nprim();
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
  auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] =
      Engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), deriv_order);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  cout << "compute_2body_fock:precision = " << precision << std::endl;
  cout << "Engine::precision = " << engines[0].precision() << std::endl;
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = obs.shell2bf();
  auto shell2atom = obs.shell2atom(atoms);

  auto lambda = [&](int thread_id) {

    auto& engine = engines[thread_id];
    const auto& buf = engine.results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto& timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    size_t shell_atoms[4];

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell
      shell_atoms[0] = shell2atom[s1];

      for (const auto& s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();
        shell_atoms[1] = shell2atom[s2];

        const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();
          shell_atoms[2] = shell2atom[s3];

          const auto Dnorm123 =
              do_schwarz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto& s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break;  // for each s3, s4 are stored in monotonically increasing
                      // order

            if ((s1234++) % nthreads != thread_id) continue;

            const auto Dnorm1234 =
                do_schwarz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwarz_screen &&
                Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();
            shell_atoms[3] = shell2atom[s4];

            const auto n1234 = n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

            // computes contribution from shell set \c idx to the operator matrix with
            // index \c op
            auto add_shellset_to_dest = [&](
                std::size_t op, std::size_t idx, int coord1, int coord2, double scale = 1.0) {
              auto& g = G[op];
              auto shset = buf[idx];
              const auto weight = scale * s1234_deg;

              for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for (auto f2 = 0; f2 != n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for (auto f3 = 0; f3 != n3; ++f3) {
                    const auto bf3 = f3 + bf3_first;
                    for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                      const auto bf4 = f4 + bf4_first;

                      const auto value = shset[f1234];
                      const auto wvalue = value * weight;

                      g(bf1, bf2) += D(bf3, bf4) * wvalue;
                      g(bf3, bf4) += D(bf1, bf2) * wvalue;
                      g(bf1, bf3) -= 0.25 * D(bf2, bf4) * wvalue;
                      g(bf2, bf4) -= 0.25 * D(bf1, bf3) * wvalue;
                      g(bf1, bf4) -= 0.25 * D(bf2, bf3) * wvalue;
                      g(bf2, bf3) -= 0.25 * D(bf1, bf4) * wvalue;
                    }
                  }
                }
              }
            };

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<Operator::coulomb, BraKet::xx_xx, deriv_order>(
                obs[s1], obs[s2], obs[s3], obs[s4]);
            if (buf[0] == nullptr)
              continue; // if all integrals screened out, skip to next quartet
            num_ints_computed += nderiv_shellset * n1234;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            switch (deriv_order) {
              case 0: {
                int coord1 = 0, coord2 = 0;
                add_shellset_to_dest(thread_id, 0, coord1, coord2);
              } break;

              case 1: {
                for (auto d = 0; d != 12; ++d) {
                  const int a = d / 3;
                  const int xyz = d % 3;

                  auto coord = shell_atoms[a] * 3 + xyz;
                  auto& g = G[thread_id * nderiv + coord];

                  int coord1 = 0, coord2 = 0;

                  add_shellset_to_dest(thread_id * nderiv + coord, d, coord1, coord2);

                }  // d \in [0,12)
              } break;

              case 2: {
// computes upper triangle index
// n2 = matrix size times 2
// i,j = (unordered) indices
#define upper_triangle_index(n2, i, j)                           \
  (std::min((i), (j))) * ((n2) - (std::min((i), (j))) - 1) / 2 + \
      (std::max((i), (j)))
                // look over shellsets in the order in which they appear
                std::size_t shellset_idx = 0;
                for (auto c1 = 0; c1 != 4; ++c1) {
                  auto a1 = shell_atoms[c1];
                  auto coord1 = 3 * a1;
                  for (auto xyz1 = 0; xyz1 != 3; ++xyz1, ++coord1) {
                    for (auto c2 = c1; c2 != 4; ++c2) {
                      auto a2 = shell_atoms[c2];
                      auto xyz2_start = (c1 == c2) ? xyz1 : 0;
                      auto coord2 = 3 * a2 + xyz2_start;
                      for (auto xyz2 = xyz2_start; xyz2 != 3;
                           ++xyz2, ++coord2) {
                        double scale =
                            (coord1 == coord2 && c1 != c2) ? 2.0 : 1.0;

                        const auto coord12 = upper_triangle_index(
                            ncoords_times_two, coord1, coord2);
                        auto op = thread_id * nderiv + coord12;
                        add_shellset_to_dest(op, shellset_idx, coord1, coord2, scale);
                        ++shellset_idx;
                      }
                    }
                  }
                }
              } break;
#undef upper_triangle_index

              default:
                assert(deriv_order <= 2 &&
                       "support for 3rd and higher derivatives of the Fock "
                       "matrix not yet implemented");
            }
          }
        }
      }
    }

  };  // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t t = 1; t != nthreads; ++t) {
    for (auto d = 0; d != nderiv; ++d) {
      G[d] += G[t * nderiv + d];
    }
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto& t : timers) {
    time_for_ints += t.read(0);
  }
  cout << "time for integrals = " << time_for_ints << std::endl;
  for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
#endif

  std::vector<Matrix> GG(nderiv);
  for (auto d = 0; d != nderiv; ++d) {
    GG[d] = 0.5 * (G[d] + G[d].transpose());
  }

  cout << "# of integrals = " << num_ints_computed << std::endl;

  // symmetrize the result and return
  return GG;
}

#endif










#if LIBINT2_DERIV_ERI_ORDER
template <unsigned deriv_order>
std::vector<VectorOfMatrices> HartreeFockClass::compute_2body_fock_uhf_deriv(const BasisSet& obs,
//std::vector<std::vector<Matrix>> HartreeFockClass::compute_2body_fock_uhf_deriv(const BasisSet& obs,
                                             const std::vector<Atom>& atoms,
                                             const Matrix& Dalpha, const Matrix& Dbeta, double precision,
                                             const Matrix& Schwarz) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  const auto nderiv_shellset =
      libint2::num_geometrical_derivatives(4, deriv_order); // # of derivs for each shell quartet
  const auto nderiv = libint2::num_geometrical_derivatives(
      atoms.size(), deriv_order);  // total # of derivs
  const auto ncoords_times_two = (atoms.size() * 3) * 2;
  using libint2::nthreads;

  std::vector<Matrix> Galpha(nthreads * nderiv, Matrix::Zero(n, n));
  std::vector<Matrix> Gbeta(nthreads * nderiv, Matrix::Zero(n, n));


  const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;

  Matrix Dalpha_shblk_norm =
      compute_shellblock_norm(obs, Dalpha);  // matrix of infty-norms of shell blocks
  Matrix Dbeta_shblk_norm =
      compute_shellblock_norm(obs, Dbeta);  // matrix of infty-norms of shell blocks


  auto fock_precision = precision;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  auto max_nprim = obs.max_nprim();
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision = std::min(fock_precision / (std::max(Dalpha_shblk_norm.maxCoeff(), Dbeta_shblk_norm.maxCoeff())),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;


  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] =
      Engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), deriv_order);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  cout << "compute_2body_fock:precision = " << precision << std::endl;
  cout << "Engine::precision = " << engines[0].precision() << std::endl;
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = obs.shell2bf();
  auto shell2atom = obs.shell2atom(atoms);

  auto lambda = [&](int thread_id) {

    auto& engine = engines[thread_id];
    const auto& buf = engine.results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto& timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    size_t shell_atoms[4];

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell
      shell_atoms[0] = shell2atom[s1];

      for (const auto& s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();
        shell_atoms[1] = shell2atom[s2];

        const auto Dalphanorm12 = do_schwarz_screen ? Dalpha_shblk_norm(s1, s2) : 0.;
        const auto Dbetanorm12 = do_schwarz_screen ? Dbeta_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();
          shell_atoms[2] = shell2atom[s3];

          const auto Dalphanorm123 =
              do_schwarz_screen
                  ? std::max(Dalpha_shblk_norm(s1, s3),
                             std::max(Dalpha_shblk_norm(s2, s3), Dalphanorm12))
                  : 0.;
          const auto Dbetanorm123 =
              do_schwarz_screen
                  ? std::max(Dbeta_shblk_norm(s1, s3),
                             std::max(Dbeta_shblk_norm(s2, s3), Dbetanorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto& s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break;  // for each s3, s4 are stored in monotonically increasing
                      // order

            if ((s1234++) % nthreads != thread_id) continue;

            const auto Dalphanorm1234 =
                do_schwarz_screen
                    ? std::max(
                          Dalpha_shblk_norm(s1, s4),
                          std::max(Dalpha_shblk_norm(s2, s4),
                                   std::max(Dalpha_shblk_norm(s3, s4), Dalphanorm123)))
                    : 0.;
            const auto Dbetanorm1234 =
                do_schwarz_screen
                    ? std::max(
                          Dbeta_shblk_norm(s1, s4),
                          std::max(Dbeta_shblk_norm(s2, s4),
                                   std::max(Dbeta_shblk_norm(s3, s4), Dbetanorm123)))
                    : 0.;

            if (do_schwarz_screen &&
                Dalphanorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();
            shell_atoms[3] = shell2atom[s4];

            const auto n1234 = n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

            // computes contribution from shell set \c idx to the operator matrix with
            // index \c op
            auto add_shellset_to_dest = [&](
                std::size_t op, std::size_t idx, int coord1, int coord2, double scale = 1.0) {
              auto& galpha = Galpha[op];
              auto& gbeta = Gbeta[op];
              auto shset = buf[idx];
              const auto weight = scale * s1234_deg;

              for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for (auto f2 = 0; f2 != n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for (auto f3 = 0; f3 != n3; ++f3) {
                    const auto bf3 = f3 + bf3_first;
                    for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                      const auto bf4 = f4 + bf4_first;

                      const auto value = shset[f1234];
                      const auto wvalue = value * weight;

                      galpha(bf1, bf2) += (Dalpha(bf3, bf4) + Dbeta(bf3, bf4)) * wvalue;
                      galpha(bf3, bf4) += (Dalpha(bf1, bf2) + Dbeta(bf1, bf2)) * wvalue;
                      galpha(bf1, bf3) -= 0.25 * Dalpha(bf2, bf4) * wvalue;
                      galpha(bf2, bf4) -= 0.25 * Dalpha(bf1, bf3) * wvalue;
                      galpha(bf1, bf4) -= 0.25 * Dalpha(bf2, bf3) * wvalue;
                      galpha(bf2, bf3) -= 0.25 * Dalpha(bf1, bf4) * wvalue;

                      gbeta(bf1, bf2) += (Dalpha(bf3, bf4) + Dbeta(bf3, bf4)) * wvalue;
                      gbeta(bf3, bf4) += (Dalpha(bf1, bf2) + Dbeta(bf1, bf2)) * wvalue;
                      gbeta(bf1, bf3) -= 0.25 * Dbeta(bf2, bf4) * wvalue;
                      gbeta(bf2, bf4) -= 0.25 * Dbeta(bf1, bf3) * wvalue;
                      gbeta(bf1, bf4) -= 0.25 * Dbeta(bf2, bf3) * wvalue;
                      gbeta(bf2, bf3) -= 0.25 * Dbeta(bf1, bf4) * wvalue;
                    }
                  }
                }
              }
            };

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<Operator::coulomb, BraKet::xx_xx, deriv_order>(
                obs[s1], obs[s2], obs[s3], obs[s4]);
            if (buf[0] == nullptr)
              continue; // if all integrals screened out, skip to next quartet
            num_ints_computed += nderiv_shellset * n1234;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            switch (deriv_order) {
              case 0: {
                int coord1 = 0, coord2 = 0;
                add_shellset_to_dest(thread_id, 0, coord1, coord2);
              } break;

              case 1: {
                for (auto d = 0; d != 12; ++d) {
                  const int a = d / 3;
                  const int xyz = d % 3;

                  auto coord = shell_atoms[a] * 3 + xyz;
                  auto& galpha = Galpha[thread_id * nderiv + coord];
                  auto& gbeta = Gbeta[thread_id * nderiv + coord];

                  int coord1 = 0, coord2 = 0;

                  add_shellset_to_dest(thread_id * nderiv + coord, d, coord1, coord2);

                }  // d \in [0,12)
              } break;

              case 2: {
// computes upper triangle index
// n2 = matrix size times 2
// i,j = (unordered) indices
#define upper_triangle_index(n2, i, j)                           \
  (std::min((i), (j))) * ((n2) - (std::min((i), (j))) - 1) / 2 + \
      (std::max((i), (j)))
                // look over shellsets in the order in which they appear
                std::size_t shellset_idx = 0;
                for (auto c1 = 0; c1 != 4; ++c1) {
                  auto a1 = shell_atoms[c1];
                  auto coord1 = 3 * a1;
                  for (auto xyz1 = 0; xyz1 != 3; ++xyz1, ++coord1) {
                    for (auto c2 = c1; c2 != 4; ++c2) {
                      auto a2 = shell_atoms[c2];
                      auto xyz2_start = (c1 == c2) ? xyz1 : 0;
                      auto coord2 = 3 * a2 + xyz2_start;
                      for (auto xyz2 = xyz2_start; xyz2 != 3;
                           ++xyz2, ++coord2) {
                        double scale =
                            (coord1 == coord2 && c1 != c2) ? 2.0 : 1.0;

                        const auto coord12 = upper_triangle_index(
                            ncoords_times_two, coord1, coord2);
                        auto op = thread_id * nderiv + coord12;
                        add_shellset_to_dest(op, shellset_idx, coord1, coord2, scale);
                        ++shellset_idx;
                      }
                    }
                  }
                }
              } break;
#undef upper_triangle_index

              default:
                assert(deriv_order <= 2 &&
                       "support for 3rd and higher derivatives of the Fock "
                       "matrix not yet implemented");
            }
          }
        }
      }
    }

  };  // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t t = 1; t != nthreads; ++t) {
    for (auto d = 0; d != nderiv; ++d) {
      Galpha[d] += Galpha[t * nderiv + d];
      Gbeta[d] += Gbeta[t * nderiv + d];
    }
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto& t : timers) {
    time_for_ints += t.read(0);
  }
  cout << "time for integrals = " << time_for_ints << std::endl;
  for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
#endif

  //std::vector<Matrix> GG(nderiv); // ORIGINAL CODE.
  std::vector<VectorOfMatrices> GG(nderiv);
  //std::vector<std::vector<Matrix>> GG(nderiv,2); // DOES NOT WORK.
  //std::vector<std::vector<Matrix>> GG(nderiv);
  for (auto d = 0; d != nderiv; ++d) {
	GG[d].resize(2);
	//cout << "START printEigenmatrixCleanFmt(Galpha["<<d<<"])\n";
	//printEigenmatrixCleanFmt(Galpha[d]);
	//cout << "END printEigenmatrixCleanFmt(Galpha["<<d<<"])\n";
	//cout << "START printEigenmatrixCleanFmt(Gbeta["<<d<<"])\n";
	//printEigenmatrixCleanFmt(Gbeta[d]);
	//cout << "END printEigenmatrixCleanFmt(Gbeta["<<d<<"])\n";

	GG[d] = {0.5 * ( Galpha[d] + Galpha[d].transpose() ) , 0.5 * ( Gbeta[d] + Gbeta[d].transpose() ) };

    //cout << "START printEigenmatrixCleanFmt(GG["<<d<<"][0])\n";
    //printEigenmatrixCleanFmt(GG[d][0]);
    //cout << "END printEigenmatrixCleanFmt(GG["<<d<<"][0])\n";
    //cout << "START printEigenmatrixCleanFmt(GG["<<d<<"][1])\n";
    //printEigenmatrixCleanFmt(GG[d][1]);
    //cout << "END printEigenmatrixCleanFmt(GG["<<d<<"][1])\n";
  }

  cout << "# of integrals = " << num_ints_computed << std::endl;

  // symmetrize the result and return
  return GG;
}

#endif











/*
#if LIBINT2_DERIV_ERI_ORDER
template <unsigned deriv_order>
std::vector<VectorOfMatrices> HartreeFockClass::compute_2body_JK_uhf_deriv(const BasisSet& obs,
//std::vector<std::vector<Matrix>> HartreeFockClass::compute_2body_fock_uhf_deriv(const BasisSet& obs,
                                             const std::vector<Atom>& atoms,
                                             const Matrix& Dalpha, const Matrix& Dbeta, double precision,
                                             const Matrix& Schwarz) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  const auto nderiv_shellset =
      libint2::num_geometrical_derivatives(4, deriv_order); // # of derivs for each shell quartet
  const auto nderiv = libint2::num_geometrical_derivatives(
      atoms.size(), deriv_order);  // total # of derivs
  const auto ncoords_times_two = (atoms.size() * 3) * 2;
  using libint2::nthreads;

  std::vector<Matrix> J(nthreads * nderiv, Matrix::Zero(n, n));
  std::vector<Matrix> Kalpha(nthreads * nderiv, Matrix::Zero(n, n));
  std::vector<Matrix> Kbeta(nthreads * nderiv, Matrix::Zero(n, n));


  const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;

  Matrix Dalpha_shblk_norm =
      compute_shellblock_norm(obs, Dalpha);  // matrix of infty-norms of shell blocks
  Matrix Dbeta_shblk_norm =
      compute_shellblock_norm(obs, Dbeta);  // matrix of infty-norms of shell blocks


  auto fock_precision = precision;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  auto max_nprim = obs.max_nprim();
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;

  auto engine_precision = std::min(fock_precision / (std::max(Dalpha_shblk_norm.maxCoeff(), Dbeta_shblk_norm.maxCoeff())),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;


  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] =
      Engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), deriv_order);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  cout << "compute_2body_fock:precision = " << precision << std::endl;
  cout << "Engine::precision = " << engines[0].precision() << std::endl;
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = obs.shell2bf();
  auto shell2atom = obs.shell2atom(atoms);

  auto lambda = [&](int thread_id) {

    auto& engine = engines[thread_id];
    const auto& buf = engine.results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto& timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    size_t shell_atoms[4];

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell
      shell_atoms[0] = shell2atom[s1];

      for (const auto& s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();
        shell_atoms[1] = shell2atom[s2];

        const auto Dalphanorm12 = do_schwarz_screen ? Dalpha_shblk_norm(s1, s2) : 0.;
        const auto Dbetanorm12 = do_schwarz_screen ? Dbeta_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();
          shell_atoms[2] = shell2atom[s3];

          const auto Dalphanorm123 =
              do_schwarz_screen
                  ? std::max(Dalpha_shblk_norm(s1, s3),
                             std::max(Dalpha_shblk_norm(s2, s3), Dalphanorm12))
                  : 0.;
          const auto Dbetanorm123 =
              do_schwarz_screen
                  ? std::max(Dbeta_shblk_norm(s1, s3),
                             std::max(Dbeta_shblk_norm(s2, s3), Dbetanorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto& s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break;  // for each s3, s4 are stored in monotonically increasing
                      // order

            if ((s1234++) % nthreads != thread_id) continue;

            const auto Dalphanorm1234 =
                do_schwarz_screen
                    ? std::max(
                          Dalpha_shblk_norm(s1, s4),
                          std::max(Dalpha_shblk_norm(s2, s4),
                                   std::max(Dalpha_shblk_norm(s3, s4), Dalphanorm123)))
                    : 0.;
            const auto Dbetanorm1234 =
                do_schwarz_screen
                    ? std::max(
                          Dbeta_shblk_norm(s1, s4),
                          std::max(Dbeta_shblk_norm(s2, s4),
                                   std::max(Dbeta_shblk_norm(s3, s4), Dbetanorm123)))
                    : 0.;

            if (do_schwarz_screen &&
                Dalphanorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) <
                    fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();
            shell_atoms[3] = shell2atom[s4];

            const auto n1234 = n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

            // computes contribution from shell set \c idx to the operator matrix with
            // index \c op
            auto add_shellset_to_dest = [&](
                std::size_t op, std::size_t idx, int coord1, int coord2, double scale = 1.0) {
              //auto& galpha = Galpha[op];
              //auto& gbeta = Gbeta[op];
              auto& j = J[op];
              auto& kalpha = Kalpha[op];
              auto& kbeta = Kbeta[op];
              auto shset = buf[idx];
              const auto weight = scale * s1234_deg;

              for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for (auto f2 = 0; f2 != n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for (auto f3 = 0; f3 != n3; ++f3) {
                    const auto bf3 = f3 + bf3_first;
                    for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                      const auto bf4 = f4 + bf4_first;

                      const auto value = shset[f1234];
                      const auto wvalue = value * weight;

                      j(bf1, bf2) += (Dalpha(bf3, bf4) + Dbeta(bf3, bf4)) * wvalue;
                      j(bf3, bf4) += (Dalpha(bf1, bf2) + Dbeta(bf1, bf2)) * wvalue;
                      kalpha(bf1, bf3) -= 0.25 * Dalpha(bf2, bf4) * wvalue;
                      kalpha(bf2, bf4) -= 0.25 * Dalpha(bf1, bf3) * wvalue;
                      kalpha(bf1, bf4) -= 0.25 * Dalpha(bf2, bf3) * wvalue;
                      kalpha(bf2, bf3) -= 0.25 * Dalpha(bf1, bf4) * wvalue;

                      //Jalpha(bf1, bf2) += (Dalpha(bf3, bf4) + Dbeta(bf3, bf4)) * wvalue;
                      //Jbeta(bf3, bf4) += (Dalpha(bf1, bf2) + Dbeta(bf1, bf2)) * wvalue;
                      kbeta(bf1, bf3) -= 0.25 * Dbeta(bf2, bf4) * wvalue;
                      kbeta(bf2, bf4) -= 0.25 * Dbeta(bf1, bf3) * wvalue;
                      kbeta(bf1, bf4) -= 0.25 * Dbeta(bf2, bf3) * wvalue;
                      kbeta(bf2, bf3) -= 0.25 * Dbeta(bf1, bf4) * wvalue;
                    }
                  }
                }
              }
            };

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<Operator::coulomb, BraKet::xx_xx, deriv_order>(
                obs[s1], obs[s2], obs[s3], obs[s4]);
            if (buf[0] == nullptr)
              continue; // if all integrals screened out, skip to next quartet
            num_ints_computed += nderiv_shellset * n1234;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            switch (deriv_order) {
              case 0: {
                int coord1 = 0, coord2 = 0;
                add_shellset_to_dest(thread_id, 0, coord1, coord2);
              } break;

              case 1: {
                for (auto d = 0; d != 12; ++d) {
                  const int a = d / 3;
                  const int xyz = d % 3;

                  auto coord = shell_atoms[a] * 3 + xyz;
                  auto& j = J[thread_id * nderiv + coord];
                  auto& kalpha = Kalpha[thread_id * nderiv + coord];
                  auto& kbeta = Kbeta[thread_id * nderiv + coord];

                  int coord1 = 0, coord2 = 0;

                  add_shellset_to_dest(thread_id * nderiv + coord, d, coord1, coord2);

                }  // d \in [0,12)
              } break;

              case 2: {
// computes upper triangle index
// n2 = matrix size times 2
// i,j = (unordered) indices
#define upper_triangle_index(n2, i, j)                           \
  (std::min((i), (j))) * ((n2) - (std::min((i), (j))) - 1) / 2 + \
      (std::max((i), (j)))
                // look over shellsets in the order in which they appear
                std::size_t shellset_idx = 0;
                for (auto c1 = 0; c1 != 4; ++c1) {
                  auto a1 = shell_atoms[c1];
                  auto coord1 = 3 * a1;
                  for (auto xyz1 = 0; xyz1 != 3; ++xyz1, ++coord1) {
                    for (auto c2 = c1; c2 != 4; ++c2) {
                      auto a2 = shell_atoms[c2];
                      auto xyz2_start = (c1 == c2) ? xyz1 : 0;
                      auto coord2 = 3 * a2 + xyz2_start;
                      for (auto xyz2 = xyz2_start; xyz2 != 3;
                           ++xyz2, ++coord2) {
                        double scale =
                            (coord1 == coord2 && c1 != c2) ? 2.0 : 1.0;

                        const auto coord12 = upper_triangle_index(
                            ncoords_times_two, coord1, coord2);
                        auto op = thread_id * nderiv + coord12;
                        add_shellset_to_dest(op, shellset_idx, coord1, coord2, scale);
                        ++shellset_idx;
                      }
                    }
                  }
                }
              } break;
#undef upper_triangle_index

              default:
                assert(deriv_order <= 2 &&
                       "support for 3rd and higher derivatives of the Fock "
                       "matrix not yet implemented");
            }
          }
        }
      }
    }

  };  // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t t = 1; t != nthreads; ++t) {
    for (auto d = 0; d != nderiv; ++d) {
      //Galpha[d] += Galpha[t * nderiv + d];
      //Gbeta[d] += Gbeta[t * nderiv + d];
      J[d] += J[t * nderiv + d];
      Kalpha[d] += Kalpha[t * nderiv + d];
      Kbeta[d] += Kbeta[t * nderiv + d];
    }
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto& t : timers) {
    time_for_ints += t.read(0);
  }
  cout << "time for integrals = " << time_for_ints << std::endl;
  for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
#endif

  //std::vector<Matrix> GG(nderiv); // ORIGINAL CODE.
  //std::vector<VectorOfMatrices> GG(nderiv);
  std::vector<VectorOfMatrices> JKvectorAB(nderiv);
  //std::vector<std::vector<Matrix>> GG(nderiv,2); // DOES NOT WORK.
  //std::vector<std::vector<Matrix>> GG(nderiv);
  for (auto d = 0; d != nderiv; ++d) {
	//GG[d].resize(2);
	  JKvectorAB[d].resize(3);

	//cout << "START printEigenmatrixCleanFmt(Galpha["<<d<<"])\n";
	//printEigenmatrixCleanFmt(Galpha[d]);
	//cout << "END printEigenmatrixCleanFmt(Galpha["<<d<<"])\n";
	//cout << "START printEigenmatrixCleanFmt(Gbeta["<<d<<"])\n";
	//printEigenmatrixCleanFmt(Gbeta[d]);
	//cout << "END printEigenmatrixCleanFmt(Gbeta["<<d<<"])\n";

	//GG[d] = {0.5 * ( Galpha[d] + Galpha[d].transpose() ) , 0.5 * ( Gbeta[d] + Gbeta[d].transpose() ) };
	JKvectorAB[d] = {0.5 * ( J[d] + J[d].transpose() ) , 0.5 * ( Kalpha[d] + Kalpha[d].transpose() ) , 0.5 * ( Kbeta[d] + Kbeta[d].transpose() ) };

    //cout << "START printEigenmatrixCleanFmt(GG["<<d<<"][0])\n";
    //printEigenmatrixCleanFmt(GG[d][0]);
    //cout << "END printEigenmatrixCleanFmt(GG["<<d<<"][0])\n";
    //cout << "START printEigenmatrixCleanFmt(GG["<<d<<"][1])\n";
    //printEigenmatrixCleanFmt(GG[d][1]);
    //cout << "END printEigenmatrixCleanFmt(GG["<<d<<"][1])\n";
  }

  cout << "# of integrals = " << num_ints_computed << std::endl;

  // symmetrize the result and return
  //return GG;
  return JKvectorAB;
}

#endif





// accumulate contributions from all threads
for (size_t i = 1; i != nthreads; ++i) {
	J[0] += J[i];
  Kalpha[0] += Kalpha[i];
  Kbeta[0] += Kbeta[i];
}

#if defined(REPORT_INTEGRAL_TIMINGS)
double time_for_ints = 0.0;
for (auto& t : timers) {
  time_for_ints += t.read(0);
}
cout << "time for integrals = " << time_for_ints << std::endl;
for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
#endif

// SYMMETRIZE RESULT.
Matrix JJ = 0.5 * (J[0] + J[0].transpose());
Matrix KKalpha = 0.5 * (Kalpha[0] + Kalpha[0].transpose());
Matrix KKbeta = 0.5 * (Kbeta[0] + Kbeta[0].transpose());

cout << "# of integrals = " << num_ints_computed << std::endl;

VectorOfMatrices JKvectorAB(3);
JKvectorAB = {JJ,KKalpha,KKbeta};
// symmetrize the result and return
return JKvectorAB;
}




END UNNECESSARY compute_2body_JK_uhf_deriv
*/





























Matrix HartreeFockClass::compute_2body_fock_general(const BasisSet& obs, const Matrix& D,
                                  const BasisSet& D_bs, bool D_is_shelldiagonal,
                                  double precision) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  const auto n_D = D_bs.nbf();
  assert(D.cols() == D.rows() && D.cols() == n_D);

  using libint2::nthreads;
  std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb,
                      std::max(obs.max_nprim(), D_bs.max_nprim()),
                      std::max(obs.max_l(), D_bs.max_l()), 0);
  engines[0].set_precision(precision);  // shellset-dependent precision control
                                        // will likely break positive
                                        // definiteness
                                        // stick with this simple recipe
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  auto shell2bf = obs.shell2bf();
  auto shell2bf_D = D_bs.shell2bf();

  auto lambda = [&](int thread_id) {

    auto& engine = engines[thread_id];
    auto& g = G[thread_id];
    const auto& buf = engine.results();

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell

      for (auto s2 = 0; s2 <= s1; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        for (auto s3 = 0; s3 < D_bs.size(); ++s3) {
          auto bf3_first = shell2bf_D[s3];
          auto n3 = D_bs[s3].size();

          auto s4_begin = D_is_shelldiagonal ? s3 : 0;
          auto s4_fence = D_is_shelldiagonal ? s3 + 1 : D_bs.size();

          for (auto s4 = s4_begin; s4 != s4_fence; ++s4, ++s1234) {
            if (s1234 % nthreads != thread_id) continue;

            auto bf4_first = shell2bf_D[s4];
            auto n4 = D_bs[s4].size();

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

            if (s3 >= s4) {
              auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
              auto s1234_deg = s12_deg * s34_deg;
              // auto s1234_deg = s12_deg;
              engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                  obs[s1], obs[s2], D_bs[s3], D_bs[s4]);
              const auto* buf_1234 = buf[0];
              if (buf_1234 != nullptr) {
                for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                  const auto bf1 = f1 + bf1_first;
                  for (auto f2 = 0; f2 != n2; ++f2) {
                    const auto bf2 = f2 + bf2_first;
                    for (auto f3 = 0; f3 != n3; ++f3) {
                      const auto bf3 = f3 + bf3_first;
                      for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                        const auto bf4 = f4 + bf4_first;

                        const auto value = buf_1234[f1234];
                        const auto value_scal_by_deg = value * s1234_deg;
                        g(bf1, bf2) += 2.0 * D(bf3, bf4) * value_scal_by_deg;
                      }
                    }
                  }
                }
              }
            }

            engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                obs[s1], D_bs[s3], obs[s2], D_bs[s4]);
            const auto* buf_1324 = buf[0];
            if (buf_1324 == nullptr)
              continue; // if all integrals screened out, skip to next quartet

            for (auto f1 = 0, f1324 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for (auto f3 = 0; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for (auto f2 = 0; f2 != n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for (auto f4 = 0; f4 != n4; ++f4, ++f1324) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf_1324[f1324];
                    const auto value_scal_by_deg = value * s12_deg;
                    g(bf1, bf2) -= D(bf3, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }

  };  // thread lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    G[0] += G[i];
  }

  // symmetrize the result and return
  return 0.5 * (G[0] + G[0].transpose());
}

#ifdef HAVE_DENSITY_FITTING

Matrix HartreeFockClass::DFFockEngine::compute_2body_fock_dfC(const Matrix& Cocc) {

  using libint2::nthreads;

  const auto n = obs.nbf();
  const auto ndf = dfbs.nbf();

  libint2::Timers<1> wall_timer;
  wall_timer.set_now_overhead(25);
  std::vector<libint2::Timers<5>> timers(nthreads);
  for(auto& timer: timers) timer.set_now_overhead(25);

  typedef btas::RangeNd<CblasRowMajor, std::array<long, 1>> Range1d;
  typedef btas::RangeNd<CblasRowMajor, std::array<long, 2>> Range2d;
  typedef btas::Tensor<double, Range1d> Tensor1d;
  typedef btas::Tensor<double, Range2d> Tensor2d;

  // using first time? compute 3-center ints and transform to inv sqrt
  // representation
  if (xyK.size() == 0) {

    wall_timer.start(0);

    const auto nshells = obs.size();
    const auto nshells_df = dfbs.size();
    const auto& unitshell = libint2::Shell::unit();

    // construct the 2-electron 3-center repulsion integrals engine
    // since the code assumes (xx|xs) braket, and Engine/libint only produces
    // (xs|xx), use 4-center engine
    std::vector<libint2::Engine> engines(nthreads);
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
                                 std::max(obs.max_nprim(), dfbs.max_nprim()),
                                 std::max(obs.max_l(), dfbs.max_l()), 0);
    engines[0].set_braket(BraKet::xs_xx);
    for (size_t i = 1; i != nthreads; ++i) {
      engines[i] = engines[0];
    }

    auto shell2bf = obs.shell2bf();
    auto shell2bf_df = dfbs.shell2bf();

    Tensor3d Zxy{ndf, n, n};

    auto lambda = [&](int thread_id) {

      auto& engine = engines[thread_id];
      auto& timer = timers[thread_id];
      const auto& results = engine.results();

      // loop over permutationally-unique set of shells
      for (auto s1 = 0l, s123 = 0l; s1 != nshells_df; ++s1) {
        auto bf1_first = shell2bf_df[s1];  // first basis function in this shell
        auto n1 = dfbs[s1].size();  // number of basis functions in this shell

        for (auto s2 = 0; s2 != nshells; ++s2) {
          auto bf2_first = shell2bf[s2];
          auto n2 = obs[s2].size();
          const auto n12 = n1 * n2;

          for (auto s3 = 0; s3 != nshells; ++s3, ++s123) {
            if (s123 % nthreads != thread_id) continue;

            auto bf3_first = shell2bf[s3];
            auto n3 = obs[s3].size();
            const auto n123 = n12 * n3;

            timer.start(0);

            engine.compute2<Operator::coulomb, BraKet::xs_xx, 0>(
                dfbs[s1], unitshell, obs[s2], obs[s3]);
            const auto* buf = results[0];
            if (buf == nullptr)
              continue;

            timer.stop(0);
            timer.start(1);

            auto lower_bound = {bf1_first, bf2_first, bf3_first};
            auto upper_bound = {bf1_first + n1, bf2_first + n2, bf3_first + n3};
            auto view = btas::make_view(
                Zxy.range().slice(lower_bound, upper_bound), Zxy.storage());
            std::copy(buf, buf + n123, view.begin());

            timer.stop(1);
          }  // s3
        }    // s2
      }      // s1

    };  // lambda

    libint2::parallel_do(lambda);

    wall_timer.stop(0);

    double ints_time = 0;
    for(const auto& timer: timers) ints_time += timer.read(0);
    cout << "time for Zxy integrals = " << ints_time << " (total from all threads)" << std::endl;
    double copy_time = 0;
    for(const auto& timer: timers) copy_time += timer.read(1);
    cout << "time for copying into BTAS = " << copy_time << " (total from all threads)"<< std::endl;
    cout << "wall time for Zxy integrals + copy = " << wall_timer.read(0) << std::endl;

    timers[0].start(2);

    Matrix V = compute_2body_2index_ints(dfbs);
    Eigen::LLT<Matrix> V_LLt(V);
    Matrix I = Matrix::Identity(ndf, ndf);
    auto L = V_LLt.matrixL();
    Matrix V_L = L;
    Matrix Linv = L.solve(I).transpose();
    // check
    //  cout << "||V - L L^t|| = " << (V - V_L * V_L.transpose()).norm() <<
    //  std::endl;
    //  cout << "||I - L L^-1^t|| = " << (I - V_L *
    //  Linv.transpose()).norm() << std::endl;
    //  cout << "||V^-1 - L^-1 L^-1^t|| = " << (V.inverse() - Linv *
    //  Linv.transpose()).norm() << std::endl;

    Tensor2d K{ndf, ndf};
    std::copy(Linv.data(), Linv.data() + ndf * ndf, K.begin());

    xyK = Tensor3d{n, n, ndf};
    btas::contract(1.0, Zxy, {1, 2, 3}, K, {1, 4}, 0.0, xyK, {2, 3, 4});
    Zxy = Tensor3d{0, 0, 0};  // release memory

    timers[0].stop(2);
    cout << "time for integrals metric tform = " << timers[0].read(2)
              << std::endl;
  }  // if (xyK.size() == 0)

  // compute exchange
  timers[0].start(3);

  const auto nocc = Cocc.cols();
  Tensor2d Co{n, nocc};
  std::copy(Cocc.data(), Cocc.data() + n * nocc, Co.begin());
  Tensor3d xiK{n, nocc, ndf};
  btas::contract(1.0, xyK, {1, 2, 3}, Co, {2, 4}, 0.0, xiK, {1, 4, 3});

  Tensor2d G{n, n};
  btas::contract(1.0, xiK, {1, 2, 3}, xiK, {4, 2, 3}, 0.0, G, {1, 4});

  timers[0].stop(3);
  cout << "time for exchange = " << timers[0].read(3) << std::endl;

  // compute Coulomb
  timers[0].start(4);

  Tensor1d Jtmp{ndf};
  btas::contract(1.0, xiK, {1, 2, 3}, Co, {1, 2}, 0.0, Jtmp, {3});
  xiK = Tensor3d{0, 0, 0};
  btas::contract(2.0, xyK, {1, 2, 3}, Jtmp, {3}, -1.0, G, {1, 2});

  timers[0].stop(4);
  cout << "time for coulomb = " << timers[0].read(4) << std::endl;

  // copy result to an Eigen::Matrix
  Matrix result(n, n);
  std::copy(G.cbegin(), G.cend(), result.data());
  return result;
}
#endif  // HAVE_DENSITY_FITTING

// should be a unit test somewhere
void HartreeFockClass::api_basic_compile_test(const BasisSet& obs) {
  using namespace libint2;
  Engine onebody_engine(
      Operator::overlap,  // will compute overlap ints
      obs.max_nprim(),    // max # of primitives in shells this engine will
                          // accept
      obs.max_l()  // max angular momentum of shells this engine will accept
      );
  auto shell2bf = obs.shell2bf();
  const auto& results = onebody_engine.results();
  for (auto s1 = 0; s1 != obs.size(); ++s1) {
    for (auto s2 = 0; s2 != obs.size(); ++s2) {
      cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
      onebody_engine.compute(obs[s1], obs[s2]);
      const auto* ints_shellset = results[0];
      cout << "done" << std::endl;

      auto bf1 = shell2bf[s1];   // first basis function in first shell
      auto n1 = obs[s1].size();  // number of basis functions in first shell
      auto bf2 = shell2bf[s2];   // first basis function in second shell
      auto n2 = obs[s2].size();  // number of basis functions in second shell

      // this iterates over integrals in the order they are packed in array
      // ints_shellset
      for (auto f1 = 0; f1 != n1; ++f1)
        for (auto f2 = 0; f2 != n2; ++f2)
          cout << "  " << bf1 + f1 << " " << bf2 + f2 << " "
                    << ints_shellset[f1 * n2 + f2] << std::endl;
    }
  }

  using libint2::Operator;

  std::vector<std::pair<double, double>> cgtg_params{
      {0.1, 0.2}, {0.3, 0.4}, {0.5, 0.6}};
  {
    auto K =
        compute_schwarz_ints<Operator::cgtg>(obs, obs, false, cgtg_params);
    cout << "cGTG Schwarz ints\n" << K << std::endl;
  }
  {
    auto K = compute_schwarz_ints<Operator::cgtg_x_coulomb>(obs, obs, false,
                                                             cgtg_params);
    cout << "cGTG/r12 Schwarz ints\n" << K << std::endl;
  }
  {
    auto K =
        compute_schwarz_ints<Operator::delcgtg2>(obs, obs, false, cgtg_params);
    cout << "||Del.cGTG||^2 Schwarz ints\n" << K << std::endl;
  }
  double attenuation_omega = 1.0;
  {
    auto K =
        compute_schwarz_ints<Operator::erfc_coulomb>(obs, obs, false, attenuation_omega);
    cout << "erfc_coulomb Schwarz ints\n" << K << std::endl;
  }
  {
    auto K =
        compute_schwarz_ints<Operator::erf_coulomb>(obs, obs, false, attenuation_omega);
    cout << "erf_coulomb Schwarz ints\n" << K << std::endl;
  }

  {  // test 2-index ints
    Engine eri4_engine(Operator::coulomb, obs.max_nprim(), obs.max_l());
    Engine eri2_engine = eri4_engine;
    eri2_engine.set_braket(BraKet::xs_xs);
    auto shell2bf = obs.shell2bf();
    const auto& results4 = eri4_engine.results();
    const auto& results2 = eri2_engine.results();
    for (auto s1 = 0; s1 != obs.size(); ++s1) {
      for (auto s2 = 0; s2 != obs.size(); ++s2) {
        eri4_engine.compute(obs[s1], Shell::unit(), obs[s2], Shell::unit());
        eri2_engine.compute(obs[s1], obs[s2]);

        auto bf1 = shell2bf[s1];   // first basis function in first shell
        auto n1 = obs[s1].size();  // number of basis functions in first shell
        auto bf2 = shell2bf[s2];   // first basis function in second shell
        auto n2 = obs[s2].size();  // number of basis functions in second shell

        const auto* buf4 = results4[0];
        const auto* buf2 = results2[0];

        // this iterates over integrals in the order they are packed in array
        // ints_shellset
        for (auto f1 = 0, f12 = 0; f1 != n1; ++f1)
          for (auto f2 = 0; f2 != n2; ++f2, ++f12)
            assert(std::abs(buf4[f12] - buf2[f12]) < 1e-12 &&
                   "2-center ints test failed");
      }
    }
  }
  {  // test 3-index ints
    Engine eri4_engine(Operator::coulomb, obs.max_nprim(), obs.max_l());
    Engine eri3_engine = eri4_engine;
    eri3_engine.set_braket(BraKet::xs_xx);
    auto shell2bf = obs.shell2bf();
    const auto& results4 = eri4_engine.results();
    const auto& results3 = eri3_engine.results();
    for (auto s1 = 0; s1 != obs.size(); ++s1) {
      for (auto s2 = 0; s2 != obs.size(); ++s2) {
        for (auto s3 = 0; s3 != obs.size(); ++s3) {
          eri4_engine.compute(obs[s1], Shell::unit(), obs[s2], obs[s3]);
          eri3_engine.compute(obs[s1], obs[s2], obs[s3]);

          auto bf1 = shell2bf[s1];   // first basis function in first shell
          auto n1 = obs[s1].size();  // number of basis functions in first shell
          auto bf2 = shell2bf[s2];   // first basis function in second shell
          auto n2 =
              obs[s2].size();       // number of basis functions in second shell
          auto bf3 = shell2bf[s3];  // first basis function in third shell
          auto n3 = obs[s3].size();  // number of basis functions in third shell

          const auto* buf4 = results4[0];
          const auto* buf3 = results3[0];

          // this iterates over integrals in the order they are packed in array
          // ints_shellset
          for (auto f1 = 0, f123 = 0; f1 != n1; ++f1)
            for (auto f2 = 0; f2 != n2; ++f2)
              for (auto f3 = 0; f3 != n3; ++f3, ++f123)
                assert(std::abs(buf4[f123] - buf3[f123]) < 1e-12 &&
                       "3-center ints test failed");
        }
      }
    }
  }

#if LIBINT2_DERIV_ERI_ORDER
  {  // test deriv 2-index ints
    Engine eri4_engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 1);
    Engine eri2_engine = eri4_engine;
    eri2_engine.set_braket(BraKet::xs_xs);
    auto shell2bf = obs.shell2bf();
    const auto& results4 = eri4_engine.results();
    const auto& results2 = eri2_engine.results();
    for (auto s1 = 0; s1 != obs.size(); ++s1) {
      for (auto s2 = 0; s2 != obs.size(); ++s2) {
        eri4_engine.compute(obs[s1], Shell::unit(), obs[s2], Shell::unit());
        eri2_engine.compute(obs[s1], obs[s2]);

        auto bf1 = shell2bf[s1];   // first basis function in first shell
        auto n1 = obs[s1].size();  // number of basis functions in first shell
        auto bf2 = shell2bf[s2];   // first basis function in second shell
        auto n2 = obs[s2].size();  // number of basis functions in second shell

        // loop over derivative shell sets
        for(auto d=0; d!=6; ++d) {
          const auto* buf4 = results4[d<3 ? d : d+3];
          const auto* buf2 = results2[d];

          // this iterates over integrals in the order they are packed in array
          // ints_shellset
          for (auto f1 = 0, f12 = 0; f1 != n1; ++f1)
            for (auto f2 = 0; f2 != n2; ++f2, ++f12)
              assert(std::abs(buf4[f12] - buf2[f12]) < 1e-12 &&
                     "deriv 2-center ints test failed");
        }

      }
    }
  }
#endif

#if LIBINT2_DERIV_ERI_ORDER > 1
  {  // test 2nd deriv 2-index ints
    Engine eri4_engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 2);
    Engine eri2_engine = eri4_engine;
    eri2_engine.set_braket(BraKet::xs_xs);
    auto shell2bf = obs.shell2bf();
    const auto& results4 = eri4_engine.results();
    const auto& results2 = eri2_engine.results();
    for (auto s1 = 0; s1 != obs.size(); ++s1) {
      for (auto s2 = 0; s2 != obs.size(); ++s2) {
        eri4_engine.compute(obs[s1], Shell::unit(), obs[s2], Shell::unit());
        eri2_engine.compute(obs[s1], obs[s2]);

        auto bf1 = shell2bf[s1];   // first basis function in first shell
        auto n1 = obs[s1].size();  // number of basis functions in first shell
        auto bf2 = shell2bf[s2];   // first basis function in second shell
        auto n2 = obs[s2].size();  // number of basis functions in second shell

        // loop over derivative shell sets
        for (auto d1 = 0, d12 = 0; d1 != 6; ++d1) {
          const auto dd1 = d1 < 3 ? d1 : d1 + 3;
          for (auto d2 = d1; d2 != 6; ++d2, ++d12) {
            const auto dd2 = d2 < 3 ? d2 : d2 + 3;
            const auto dd12 = dd1 * (24 - dd1 - 1) / 2 + dd2;
            const auto* buf4 = results4[dd12];
            const auto* buf2 = results2[d12];

            // this iterates over integrals in the order they are packed in
            // array
            // ints_shellset
            for (auto f1 = 0, f12 = 0; f1 != n1; ++f1)
              for (auto f2 = 0; f2 != n2; ++f2, ++f12)
                assert(std::abs(buf4[f12] - buf2[f12]) < 1e-12 &&
                       "2nd deriv 2-center ints test failed");
          }
        }
      }
    }
  }
#endif

}

















//****************************************************************************************************************************************
//************************************** DFT PART FROM ERKALE ****************************************************************************
//****************************************************************************************************************************************

SCF HartreeFockClass::initialize_Erkale_SCF(string basename, ErkaleBasisSet& basisp, Settings& setp, unsigned long int Nbf, int Nel, int mult, Matrix& S, Matrix& T, Matrix& Vnuc, Matrix& Hcore )
{
	  cout << "\nSTART initialize_Erkale_SCF\n";

	  // Parse settings
	  Settings set;
	  set.add_scf_settings();
	  set.add_string("SaveChk","File to use as checkpoint","erkale.chk");
	  set.add_string("LoadChk","File to load old results from","");
	  set.add_bool("ForcePol","Force polarized calculation",false);
	  //set.parse(std::string(argv[1]),true);
	  //set.parse("erkale.inp",true);
	  set.parse(basename+".inp");

	  set.add_dft_settings(); // ADDED BY C. SALGADO TO AVOID ERROR LIKE "The boolean setting DFTLobatto was not found!".

      set.print();

      // CHK part added here by C.SALGADO
      //std::string chkf(set.get_string("LoadChk"));
      //Checkpoint chk(chkf,false);
      std::string savename=set.get_string("SaveChk");
      //Checkpoint chk(savename,false);


      //chkpt.write(basis);

	  // Basis set
	  ErkaleBasisSet basis;
      //ErkaleBasisSet basis;
	  std::string basfile(set.get_string("Basis"));
	  //if(stricmp(basfile,"Read")==0) {
	  //if(strcasecmp(basfile,"Read")==0) {
	  /*
	  if(strcasecmp(basfile.c_str(),"Read")==0) {
	    // Get checkpoint file
	    std::string chkf(set.get_string("LoadChk"));
	    if(!chkf.size())
	      throw std::runtime_error("Must specify LoadChk for Basis Read\n");
	    if(!file_exists(chkf))
	      throw std::runtime_error("Can't find LoadChk!\n");
	    //Checkpoint chk(chkf,false);
	    chk.read(basis);

	    printf("Basis set read in from checkpoint.\n");

	  } else {
	  */
	    // Read in atoms.
	    std::vector<atom_t> atoms;
	    std::string atomfile=set.get_string("System");
	    if(file_exists(atomfile))
	      atoms=load_xyz(atomfile,!set.get_bool("InputBohr"));
	    else {
	      // Check if a directory has been set
	      char * libloc=getenv("ERKALE_SYSDIR");
	      if(libloc) {
	    	  //std::string filename=std::string(libloc)+"/"+atomfile;
	    	  std::string filename=atomfile;
	    	  if(file_exists(filename))
	    		  atoms=load_xyz(filename,!set.get_bool("InputBohr"));
	    	  else
	    		  throw std::runtime_error("Unable to open xyz input file!\n");
	      } else
	    	  throw std::runtime_error("Unable to open xyz input file!\n");
	    }

	    cout << "\nErkaleBasisSetLibrary baslib;\n";
	    // Read in basis set
	    ErkaleBasisSetLibrary baslib;
	    baslib.load_basis(basfile);

	    cout << "\nBEFORE construct_basis(*basis,atoms,baslib,set);\n";
	    // Construct basis set
	    construct_basis(basis,atoms,baslib,set); // NECESARIO PONER EL PUNTERO EN BASIS PARA QUE LO PILLE BIEN. VER MAIN.CPP DE ERKALE.

	    cout << "\n******************************************************\n";
	    cout << "\nPRINT BASIS in initialize_Erkale_SCF\n";
	    cout << "\n******************************************************\n";
	    basis.print(true);
	    cout << "\n******************************************************\n";
	    cout << "\nEND PRINT BASIS in initialize_Erkale_SCF\n";
	    cout << "\n******************************************************\n";
	  //}

	  //ofstream chkfile;
	  //chkfile.open(savename);

	  // Write checkpoint.
	  Checkpoint chk(savename,true); // NECESITAMOS EL TRUE CUANDO NO CARGAMOS LA BASIS DESDE UN CHK PREEXISTENTE. LÍNEA 1873 DE scf-base.cpp en Erkale.
	  chk.write(basis); // NECESARIO PONER EL PUNTERO EN BASIS PARA QUE LO PILLE BIEN. VER MAIN.CPP DE ERKALE.




	  // THIS NEEDS TO BE BEFORE INITIALIZATION OF SCF ErkaleSCF(basis,set,chk);
	  /// Overlap matrix
	  cout << "START printEigenmatrixCleanFmt(S)\n";
	  printEigenmatrixCleanFmt(S);
	  cout << "END printEigenmatrixCleanFmt(S)\n";
	  arma::mat SS = example_cast_arma(S);
	  //ErkaleSCF.set_S(example_cast_arma(S));
	  basis.set_S(SS);
	  SS.print("SS");
	  printf("\n");
	  /// Kinetic energy matrix
	  cout << "START printEigenmatrixCleanFmt(T)\n";
	  //printEigenmatrixCleanFmt(T);
	  cout << "END printEigenmatrixCleanFmt(T)\n";
	  arma::mat TT = example_cast_arma(T);
	  //ErkaleSCF.set_T(example_cast_arma(T));
	  basis.set_T(TT);
	  TT.print("TT");
	  printf("\n");
	  /// Nuclear attraction matrix
	  cout << "START printEigenmatrixCleanFmt(Vnuc)\n";
	  //printEigenmatrixCleanFmt(Vnuc);
	  cout << "END printEigenmatrixCleanFmt(Vnuc)\n";
	  arma::mat VVnuc = example_cast_arma(Vnuc);
	  //ErkaleSCF.set_Vnuc(example_cast_arma(Vnuc));
	  basis.set_Vnuc(VVnuc);
	  VVnuc.print("VVnuc");
	  printf("\n");
	  /// Core Hamiltonian
	  cout << "START printEigenmatrixCleanFmt(Hcore)\n";
	  //printEigenmatrixCleanFmt(Hcore);
	  cout << "END printEigenmatrixCleanFmt(Hcore)\n";
	  arma::mat HHcore = example_cast_arma(Hcore);
	  //ErkaleSCF.set_Hcore(example_cast_arma(Hcore));
	  basis.set_Hcore(HHcore);
	  SS.print("HHcore");
	  printf("\n");

	  //SCF(basis,set,chk) ErkaleSCF; // DOES NOT WORK.
	  SCF ErkaleSCF(basis,set,chk); // NECESARIO PONER EL PUNTERO EN BASIS PARA QUE LO PILLE BIEN. VER MAIN.CPP DE ERKALE. PARECE QUE YA NO. FIJARSE EN Matrx& NECESARIO PARA QUE MODIFICE EL CONTENIDO EN ESAS DIRECCIONES DE MEMORIA.



	  // /// Basis set to use (needed for DFT grid operation)
	  // const ErkaleBasisSet * basisp;
	  // /// Density fitting basis
	  // ErkaleBasisSet dfitbas;
	  // /// Checkpoint file
	  // Checkpoint * chkptp;

	  // /// Basis orthogonalizing matrix
	  // arma::mat Sinvh;

	  // /// Amount of basis functions
	  // size_t Nbf;


	  // /// Total number of electrons
	  // int Nel;


	  // /// Multiplicity
	  // int mult;

	  //basisp = &basis;
	  basisp = basis;
	  setp = set;

	  cout << "\n******************************************************\n";
	  cout << "\nPRINT BASISP\n";
	  cout << "\n******************************************************\n";
	  //basisp->print(true);
	  basisp.print(true);
	  cout << "\n******************************************************\n";
	  cout << "\nEND PRINT BASISP\n";
	  cout << "\n******************************************************\n";

	  cout << "END initialize_Erkale_SCF\n";

	  return ErkaleSCF;

}

//void HartreeFockClass::wrap_Fock_UDFT(string basename, const double Ea, const double Eb, const int naocc, const int nbocc, const Matrix& Calpha, const Matrix& Cbeta, const Matrix& Calpha_occ, const Matrix& Cbeta_occ, const Matrix& Halpha, const Matrix& Hbeta, const Matrix& Dalpha, const Matrix& Dbeta, const Matrix& J, const Matrix& Kalpha, const Matrix& Kbeta, const Matrix& XCalpha, const Matrix& XCbeta)
void HartreeFockClass::wrap_Fock_UDFT(SCF ErkaleSCF, ErkaleBasisSet& basis, Settings& set, string basename, unsigned long int Nbf, std::vector<double> Ea, std::vector<double> Eb, int naocc, int nbocc, Matrix& S, Matrix& Calpha, Matrix& Cbeta, Matrix& Calpha_occ, Matrix& Cbeta_occ, Matrix& Halpha, Matrix& Hbeta, Matrix& Dalpha, Matrix& Dbeta, Matrix& J, Matrix& Kalpha, Matrix& Kbeta, Matrix& XCalpha, Matrix& XCbeta)
{

	cout << "START wrap_Fock_UDFT\n";

	//Timer t;

	//ErkaleBasisSet basis = ErkaleSCF.basisp();
	//ErkaleBasisSet basis;
	//basis = ErkaleSCF.get_basis();
	//new SCF::uscf_t sol;

	arma::mat Sarma = example_cast_arma(S);

	uscf_t sol;
	  /// Orbitals
	//  arma::mat Ca, Cb;
	sol.Ca = example_cast_arma(Calpha);
	sol.Cb = example_cast_arma(Cbeta);
	arma::mat Ca_occ = example_cast_arma(Calpha_occ);
	arma::mat Cb_occ = example_cast_arma(Cbeta_occ);
	  /// Orbital energies
	//  arma::vec Ea, Eb;

	cout << "START printEigenmatrixCleanFmt(Halpha)\n";
	printEigenmatrixCleanFmt(Halpha);
	cout << "END printEigenmatrixCleanFmt(Halpha)\n";
	  /// Fock operators
	sol.Ha = example_cast_arma(Halpha);
	sol.Hb = example_cast_arma(Hbeta);
	printf("arma Ha = \n");
	sol.Ha.print("Ha");
	printf("\n");

	/// Density matrices
	sol.Pa = example_cast_arma(Dalpha);
	sol.Pb = example_cast_arma(Dbeta);
	sol.P = sol.Pa + sol.Pb;
	cout << "START printEigenmatrixCleanFmt(Dalpha)\n";
	printEigenmatrixCleanFmt(Dalpha);
	cout << "END printEigenmatrixCleanFmt(Dalpha)\n";
	printf("arma Pa = \n");
	sol.Pa.print("Pa");
	printf("\n");
	cout << "START printEigenmatrixCleanFmt(Dbeta)\n";
	printEigenmatrixCleanFmt(Dbeta);
	cout << "END printEigenmatrixCleanFmt(Dbeta)\n";
	printf("arma Pb = \n");
	sol.Pb.print("Pb");
	printf("\n");

	  /// Coulomb operator
	sol.J = example_cast_arma(J);
	  /// Exchange operators
	sol.Ka = example_cast_arma(Kalpha);
	sol.Kb = example_cast_arma(Kbeta);
	  /// KS-XC matrix
	sol.XCa = example_cast_arma(XCalpha);
	sol.XCb = example_cast_arma(XCbeta);

	sol.Ea = arma::vec(Ea);
	sol.Eb = arma::vec(Eb);



	bool verbose = true;

	/*
	// SOLVE ERROR WITH THE VALUES OF Ea AND Eb WITH OUTPUT
	// Computing DFT exchange and correlation ... done (0.00 s)
	// Numerically integrated density is 0.00000 (-100.0000 %).
	// Warning - numerically integrated density seems inaccurate.
	bool verbose = false;
	//double trrhmins = 8.673617e-18;
    double trrhmins = 1e-5;

	arma::mat Canew;
	arma::vec Eanew;
	//TRRH_update(sol.Ha,sol.Ca,S,Canew,Eanew,nocca,verbose,trrhmins);
	TRRH_update(sol.Ha,sol.Ca,Sarma,Canew,Eanew,naocc,verbose,trrhmins);

	arma::mat Cbnew;
	arma::vec Ebnew;
	//TRRH_update(sol.Hb,sol.Cb,S,Cbnew,Ebnew,noccb,verbose,trrhmins);
	TRRH_update(sol.Hb,sol.Cb,Sarma,Cbnew,Ebnew,nbocc,verbose,trrhmins);


	sol.Ca=Canew;
	sol.Cb=Cbnew;

	sol.Ea=Eanew;
	sol.Eb=Ebnew;
	// END SOLVE ERROR WITH THE VALUES OF Ea AND Eb WITH OUTPUT
	*/

	  /// Imaginary exchange
	  //Matrix Ka_im, Kb_im;

	  /// Complex orbitals (for SIC)
	  //Matrix cCa, cCb;
	  /// Imaginary part of complex-CMO density matrix (for complex exchange contribution)
	 // Matrix Pa_im, Pb_im;

	  /// Energy information
	  //energy_t en;
	//sol.en.Eel = Ea + Eb;
	sol.en.Eel = 0.0;


	// FROM MAIN.CPP
	/*
	  // Parse settings
	  Settings set;
	  set.add_scf_settings();
	  set.add_string("SaveChk","File to use as checkpoint","erkale.chk");
	  set.add_string("LoadChk","File to load old results from","");
	  set.add_bool("ForcePol","Force polarized calculation",false);
	  //set.parse(std::string(argv[1]),true);
	  //set.parse("erkale.inp",true);
	  set.parse(basename+".inp");

	  set.add_dft_settings(); // ADDED BY C. SALGADO TO AVOID ERROR LIKE "The boolean setting DFTLobatto was not found!".

	  set.print();
	  */

/*
	  // Basis set
	  ErkaleBasisSet basis;
	  std::string basfile(set.get_string("Basis"));
	  //if(stricmp(basfile,"Read")==0) {
	  //if(strcasecmp(basfile,"Read")==0) {
	  if(strcasecmp(basfile.c_str(),"Read")==0) {
	    // Get checkpoint file
	    std::string chkf(set.get_string("LoadChk"));
	    if(!chkf.size())
	      throw std::runtime_error("Must specify LoadChk for Basis Read\n");
	    if(!file_exists(chkf))
	      throw std::runtime_error("Can't find LoadChk!\n");
	    Checkpoint chk(chkf,false);
	    chk.read(basis);

	    printf("Basis set read in from checkpoint.\n");

	  } else {
	    // Read in atoms.
	    std::vector<atom_t> atoms;
	    std::string atomfile=set.get_string("System");
	    if(file_exists(atomfile))
	      atoms=load_xyz(atomfile,!set.get_bool("InputBohr"));
	    else {
	      // Check if a directory has been set
	      char * libloc=getenv("ERKALE_SYSDIR");
	      if(libloc) {
		std::string filename=std::string(libloc)+"/"+atomfile;
		if(file_exists(filename))
		  atoms=load_xyz(filename,!set.get_bool("InputBohr"));
		else
		  throw std::runtime_error("Unable to open xyz input file!\n");
	      } else
		throw std::runtime_error("Unable to open xyz input file!\n");
	    }

	    // Read in basis set
	    ErkaleBasisSetLibrary baslib;
	    baslib.load_basis(basfile);

	    // Construct basis set
	    construct_basis(basis,atoms,baslib,set);
	  }
*/

	//verbose = true;
	//dft_t dft;
	//dft=parse_dft(set,false);
	dft_t dft0;
	dft0=parse_dft(set,false);
	//dft0=parse_dft(ErkaleSCF.set,false);
	dft_t dft(dft0);

	parse_xc_func(dft.x_func,dft.c_func,set.get_string("Method"));
	//parse_xc_func(dft.x_func,dft.c_func,"UDFT");
	//parse_xc_func(dft.x_func,dft.c_func,"LDA-LDA");
	//dft=parse_dft(set,false);


    cout << "\n******************************************************\n";
    cout << "\nPRINT BASIS in wrap_Fock_UDFT\n";
    cout << "\n******************************************************\n";
    basis.print(true);
    cout << "\n******************************************************\n";
    cout << "\nEND PRINT BASIS in wrap_Fock_UDFT\n";
    cout << "\n******************************************************\n";
	//DFTGrid grid(basis,verbose,dft.lobatto);
	//DFTGrid nlgrid(basis,verbose,dft.lobatto);
	DFTGrid grid(&basis,verbose,dft.lobatto);
	DFTGrid nlgrid(&basis,verbose,dft.lobatto);
    //DFTGrid *grid = new DFTGrid(&basis,verbose,dft.lobatto);
    //DFTGrid *nlgrid = new DFTGrid(&basis,verbose,dft.lobatto);

	bool strictint = ErkaleSCF.get_strictint();
	  /*
	  if(verbose && maxiter>0) {
    fprintf(stderr,"Running ");
# 238 "/home/carlos/eclipse-workspace/2018-11-12-erkale/erkale-master/src/scf-solvers.cpp.in"
    fprintf(stderr,"unrestricted ");
    if(dft.c_func>0) {
      fprintf(stderr,"%s-%s ",get_keyword(dft.x_func).c_str(),get_keyword(dft.c_func).c_str());
    } else
      fprintf(stderr,"%s ",get_keyword(dft.x_func).c_str());
    fprintf(stderr,"calculation");
    if(densityfit) {
      fprintf(stderr," with density fitting");
    }
# 272 "/home/carlos/eclipse-workspace/2018-11-12-erkale/erkale-master/src/scf-solvers.cpp.in"
    fprintf(stderr,".\n%4s %16s %10s %9s %9s %9s %10s\n","iter","E","dE","RMS dens","MAX dens","DIIS err","titer (s)");
  }
  fflush(stdout);
  */


  if(dft.x_func>0 || dft.c_func>0) {
    if(dft.adaptive) {
      grid.construct(sol.Pa,sol.Pb,dft.gridtol,dft.x_func,dft.c_func);
    } else {
      grid.construct(dft.nrad,dft.lmax,dft.x_func,dft.c_func,strictint);
      if(dft.nl)
        nlgrid.construct(dft.nlnrad,dft.nllmax,true,false,false,strictint,true);
    }
    if(verbose) {
      fflush(stdout);
      //fprintf(stderr,"%-65s %10.3f\n","    DFT grid formation",t.get());
    }
  }


	cout << "\n******************************************************\n";
	cout << "\nPRINT GRID in wrap_Fock_UDFT\n";
	cout << "\n******************************************************\n";
	grid.print_grid(set.get_string("Method"));
	cout << "\n******************************************************\n";
	cout << "\nEND PRINT GRID in wrap_Fock_UDFT\n";
	cout << "\n******************************************************\n";

	//std::vector<double> vector2(length, 0.0);
	std::vector<double> occa(naocc,1.0), occb(nbocc,1.0);

	verbose = false;
	cout << "\nBEFORE ENTERING Fock_UDFT\n";
	//SCF::Fock_UDFT(sol,occa,occb,dft,grid,nlgrid);
	ErkaleSCF.Fock_UDFT(sol,occa,occb,dft,grid,nlgrid);
	cout << "\nAFTER EXITING Fock_UDFT\n";

	verbose = true;


	//sol.Ca.print(cout);
	  /// Orbitals
	//  arma::mat Ca, Cb;
	Calpha = example_cast_eigen(sol.Ca);
	Cbeta = example_cast_eigen(sol.Cb);
	Calpha_occ = example_cast_eigen(Ca_occ);
	Cbeta_occ = example_cast_eigen(Cb_occ);

	 /// Orbital energies
	//  arma::vec Ea, Eb;

	  /// Fock operators
	Halpha = example_cast_eigen(sol.Ha);
	Hbeta = example_cast_eigen(sol.Hb);
	  /// Density matrices
	Dalpha = example_cast_eigen(sol.Pa);
	Dbeta = example_cast_eigen(sol.Pb);
	//P = sol.Pa + sol.Pb;

	  /// Coulomb operator
	J = example_cast_eigen(sol.J);
	  /// Exchange operators
	Kalpha = example_cast_eigen(sol.Ka);
	Kbeta = example_cast_eigen(sol.Kb);
	  /// KS-XC matrix
	XCalpha = example_cast_eigen(sol.XCa);
	XCbeta = example_cast_eigen(sol.XCb);

	  /// Imaginary exchange
	  //Matrix Ka_im, Kb_im;

	  /// Complex orbitals (for SIC)
	  //Matrix cCa, cCb;
	  /// Imaginary part of complex-CMO density matrix (for complex exchange contribution)
	  // Matrix Pa_im, Pb_im;

	  /// Energy information
	  //energy_t en;
	//sol.en = Ea + Eb;
	//Ea = sol.en/2;
	//Eb = sol.en/2;
	cout << "END wrap_Fock_UDFT\n";
}




//arma::vec SCF::force_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol)
//void HartreeFockClass::wrap_Fock_UDFT(string basename, const double Ea, const double Eb, const int naocc, const int nbocc, const Matrix& Calpha, const Matrix& Cbeta, const Matrix& Calpha_occ, const Matrix& Cbeta_occ, const Matrix& Halpha, const Matrix& Hbeta, const Matrix& Dalpha, const Matrix& Dbeta, const Matrix& J, const Matrix& Kalpha, const Matrix& Kbeta, const Matrix& XCalpha, const Matrix& XCbeta)
std::vector<std::vector<Matrix> > HartreeFockClass::wrap_force_UDFT(SCF ErkaleSCF, ErkaleBasisSet& basis, Settings& set, string basename, unsigned long int Nbf, std::vector<double> Ea, std::vector<double> Eb, int naocc, int nbocc, Matrix& S, Matrix& Calpha, Matrix& Cbeta, Matrix& Calpha_occ, Matrix& Cbeta_occ, Matrix& Halpha, Matrix& Hbeta, Matrix& Dalpha, Matrix& Dbeta, Matrix& J, Matrix& Kalpha, Matrix& Kbeta, Matrix& XCalpha, Matrix& XCbeta)
{

	cout << "START wrap_force_UDFT\n";

	arma::mat Sarma = example_cast_arma(S);

	//ErkaleBasisSet basis = ErkaleSCF.basisp();
	//ErkaleBasisSet basis;
	//basis = ErkaleSCF.get_basis();
	//new SCF::uscf_t sol;
	uscf_t sol;
	  /// Orbitals
	//  arma::mat Ca, Cb;
	sol.Ca = example_cast_arma(Calpha);
	sol.Cb = example_cast_arma(Cbeta);
	arma::mat Ca_occ = example_cast_arma(Calpha_occ);
	arma::mat Cb_occ = example_cast_arma(Cbeta_occ);
	  /// Orbital energies
	//  arma::vec Ea, Eb;

	cout << "START printEigenmatrixCleanFmt(Halpha)\n";
	printEigenmatrixCleanFmt(Halpha);
	cout << "END printEigenmatrixCleanFmt(Halpha)\n";
	  /// Fock operators
	sol.Ha = example_cast_arma(Halpha);
	sol.Hb = example_cast_arma(Hbeta);
	printf("arma Ha = \n");
	sol.Ha.print("Ha");
	printf("\n");
	  /// Density matrices
	sol.Pa = example_cast_arma(Dalpha);
	sol.Pb = example_cast_arma(Dbeta);
	sol.P = sol.Pa + sol.Pb;

	  /// Coulomb operator
	sol.J = example_cast_arma(J);
	  /// Exchange operators
	sol.Ka = example_cast_arma(Kalpha);
	sol.Kb = example_cast_arma(Kbeta);
	  /// KS-XC matrix
	sol.XCa = example_cast_arma(XCalpha);
	sol.XCb = example_cast_arma(XCbeta);

	sol.Ea = arma::vec(Ea);
	sol.Eb = arma::vec(Eb);

	bool verbose = true;
	/*
	// SOLVE ERROR WITH THE VALUES OF Ea AND Eb WITH OUTPUT
	// Computing DFT exchange and correlation ... done (0.00 s)
	// Numerically integrated density is 0.00000 (-100.0000 %).
	// Warning - numerically integrated density seems inaccurate.
	bool verbose = false;
	//double trrhmins = 8.673617e-18;
	double trrhmins = 1e-5;

	arma::mat Canew;
	arma::vec Eanew;
	//TRRH_update(sol.Ha,sol.Ca,S,Canew,Eanew,nocca,verbose,trrhmins);
	TRRH_update(sol.Ha,sol.Ca,Sarma,Canew,Eanew,naocc,verbose,trrhmins);

	arma::mat Cbnew;
	arma::vec Ebnew;
	//TRRH_update(sol.Hb,sol.Cb,S,Cbnew,Ebnew,noccb,verbose,trrhmins);
	TRRH_update(sol.Hb,sol.Cb,Sarma,Cbnew,Ebnew,nbocc,verbose,trrhmins);


	sol.Ca=Canew;
	sol.Cb=Cbnew;

	sol.Ea=Eanew;
	sol.Eb=Ebnew;
	// END SOLVE ERROR WITH THE VALUES OF Ea AND Eb WITH OUTPUT
    */


	  /// Imaginary exchange
	  //Matrix Ka_im, Kb_im;

	  /// Complex orbitals (for SIC)
	  //Matrix cCa, cCb;
	  /// Imaginary part of complex-CMO density matrix (for complex exchange contribution)
	 // Matrix Pa_im, Pb_im;

	  /// Energy information
	  //energy_t en;
	//sol.en.Eel = Ea + Eb;
	sol.en.Eel = 0.0;


	// FROM MAIN.CPP
	  // Parse settings
	  /*
	  Settings set;
	  set.add_scf_settings();
	  set.add_string("SaveChk","File to use as checkpoint","erkale.chk");
	  set.add_string("LoadChk","File to load old results from","");
	  set.add_bool("ForcePol","Force polarized calculation",false);
	  //set.parse(std::string(argv[1]),true);
	  //set.parse("erkale.inp",true);
	  set.parse(basename+".inp");

	  set.add_dft_settings(); // ADDED BY C. SALGADO TO AVOID ERROR LIKE "The boolean setting DFTLobatto was not found!".

	  set.print();
	  */

/*
	  // Basis set
	  ErkaleBasisSet basis;
	  std::string basfile(set.get_string("Basis"));
	  //if(stricmp(basfile,"Read")==0) {
	  //if(strcasecmp(basfile,"Read")==0) {
	  if(strcasecmp(basfile.c_str(),"Read")==0) {
	    // Get checkpoint file
	    std::string chkf(set.get_string("LoadChk"));
	    if(!chkf.size())
	      throw std::runtime_error("Must specify LoadChk for Basis Read\n");
	    if(!file_exists(chkf))
	      throw std::runtime_error("Can't find LoadChk!\n");
	    Checkpoint chk(chkf,false);
	    chk.read(basis);

	    printf("Basis set read in from checkpoint.\n");

	  } else {
	    // Read in atoms.
	    std::vector<atom_t> atoms;
	    std::string atomfile=set.get_string("System");
	    if(file_exists(atomfile))
	      atoms=load_xyz(atomfile,!set.get_bool("InputBohr"));
	    else {
	      // Check if a directory has been set
	      char * libloc=getenv("ERKALE_SYSDIR");
	      if(libloc) {
		std::string filename=std::string(libloc)+"/"+atomfile;
		if(file_exists(filename))
		  atoms=load_xyz(filename,!set.get_bool("InputBohr"));
		else
		  throw std::runtime_error("Unable to open xyz input file!\n");
	      } else
		throw std::runtime_error("Unable to open xyz input file!\n");
	    }

	    // Read in basis set
	    ErkaleBasisSetLibrary baslib;
	    baslib.load_basis(basfile);

	    // Construct basis set
	    construct_basis(basis,atoms,baslib,set);
	  }
*/

	//verbose = true;
		//dft_t dft;
		//dft=parse_dft(set,false);
		dft_t dft0;
		dft0=parse_dft(set,false);
		//dft0=parse_dft(ErkaleSCF.set,false);
		dft_t dft(dft0);

		parse_xc_func(dft.x_func,dft.c_func,set.get_string("Method"));
		//parse_xc_func(dft.x_func,dft.c_func,"UDFT");
		//parse_xc_func(dft.x_func,dft.c_func,"LDA-LDA");
		//dft=parse_dft(set,false);


	    cout << "\n******************************************************\n";
	    cout << "\nPRINT BASIS in wrap_Fock_UDFT\n";
	    cout << "\n******************************************************\n";
	    basis.print(true);
	    cout << "\n******************************************************\n";
	    cout << "\nEND PRINT BASIS in wrap_Fock_UDFT\n";
	    cout << "\n******************************************************\n";
		//DFTGrid grid(basis,verbose,dft.lobatto);
		//DFTGrid nlgrid(basis,verbose,dft.lobatto);
		DFTGrid grid(&basis,verbose,dft.lobatto);
		DFTGrid nlgrid(&basis,verbose,dft.lobatto);
	    //DFTGrid *grid = new DFTGrid(&basis,verbose,dft.lobatto);
	    //DFTGrid *nlgrid = new DFTGrid(&basis,verbose,dft.lobatto);

		bool strictint = ErkaleSCF.get_strictint();
		  /*
		  if(verbose && maxiter>0) {
	    fprintf(stderr,"Running ");
	# 238 "/home/carlos/eclipse-workspace/2018-11-12-erkale/erkale-master/src/scf-solvers.cpp.in"
	    fprintf(stderr,"unrestricted ");
	    if(dft.c_func>0) {
	      fprintf(stderr,"%s-%s ",get_keyword(dft.x_func).c_str(),get_keyword(dft.c_func).c_str());
	    } else
	      fprintf(stderr,"%s ",get_keyword(dft.x_func).c_str());
	    fprintf(stderr,"calculation");
	    if(densityfit) {
	      fprintf(stderr," with density fitting");
	    }
	# 272 "/home/carlos/eclipse-workspace/2018-11-12-erkale/erkale-master/src/scf-solvers.cpp.in"
	    fprintf(stderr,".\n%4s %16s %10s %9s %9s %9s %10s\n","iter","E","dE","RMS dens","MAX dens","DIIS err","titer (s)");
	  }
	  fflush(stdout);
	  */


	  if(dft.x_func>0 || dft.c_func>0) {
	    if(dft.adaptive) {
	      grid.construct(sol.Pa,sol.Pb,dft.gridtol,dft.x_func,dft.c_func);
	    } else {
	      grid.construct(dft.nrad,dft.lmax,dft.x_func,dft.c_func,strictint);
	      if(dft.nl)
	        nlgrid.construct(dft.nlnrad,dft.nllmax,true,false,false,strictint,true);
	    }
	    if(verbose) {
	      fflush(stdout);
	      //fprintf(stderr,"%-65s %10.3f\n","    DFT grid formation",t.get());
	    }
	  }


		cout << "\n******************************************************\n";
		cout << "\nPRINT GRID in wrap_force_UDFT\n";
		cout << "\n******************************************************\n";
		grid.print_grid(set.get_string("Method"));
		cout << "\n******************************************************\n";
		cout << "\nEND PRINT GRID in wrap_force_UDFT\n";
		cout << "\n******************************************************\n";

		//std::vector<double> vector2(length, 0.0);
		std::vector<double> occa(naocc,1.0), occb(nbocc,1.0);


	std::vector< std::vector<arma::mat> > fmat;
	//fmat.resize(3*basis->get_Nnuc());
	fmat.resize(3*basis.get_Nnuc());
	//for(size_t i=0;i<3*basis->get_Nnuc();i++) {

	for(size_t i=0;i<3*basis.get_Nnuc();i++) {
		fmat[i].resize(2);
		fmat[i][0].resize(Calpha.rows(),Calpha.cols());
		fmat[i][1].resize(Calpha.rows(),Calpha.cols());
	}


	cout << "\nBEFORE ENTERING set_basis\n";
	ErkaleSCF.set_basis(&basis);
	cout << "\nAFTER EXITING set_basis\n";

	verbose = false;

	cout << "\nBEFORE ENTERING force_UDFT_ANT\n";
	//SCF::Fock_UDFT(sol,occa,occb,dft,grid,nlgrid);
	//ErkaleSCF.Fock_UDFT(sol,occa,occb,dft,grid,nlgrid);
	arma::vec ErkaleForceVec = ErkaleSCF.force_UDFT_ANT(sol,occa,occb,dft,grid,nlgrid,set.get_double("DFTFinalTol"),fmat);//, double tol) //DFTInitialTol
	cout << "\nAFTER EXITING force_UDFT_ANT\n";

	//Matrix ErkaleForceMat arma::reshape(ErkaleForceVec,3,basis.get_Nnuc()); // UTIL SI QUISIERAMOS DEVOLVER UNA MATRIZ, PERO DEVOLVEMOS UN VECTOR DE 3*basis.get_Nnuc() COORDENADAS.
	//Matrix F2Erk = example_cast_eigen(ErkaleForceMat);

	Matrix F2Erk = example_cast_eigen(ErkaleForceVec);

	std::vector<std::vector<Matrix> > G1Erk;
	//G1Erk.resize(3*basis->get_Nnuc());
	G1Erk.resize(3*basis.get_Nnuc());
	//for(size_t i=0;i<3*basis->get_Nnuc();i++) {
	for(size_t i=0;i<3*basis.get_Nnuc();i++) {
		G1Erk[i].resize(2);
		G1Erk[i][0].resize(Calpha.rows(),Calpha.cols());
		G1Erk[i][1].resize(Calpha.rows(),Calpha.cols());
		G1Erk[i][0] = example_cast_eigen(fmat[i][0]);
		G1Erk[i][1] = example_cast_eigen(fmat[i][1]);
	}

	verbose = true;

	  /// Orbitals
	//  arma::mat Ca, Cb;
	Calpha = example_cast_eigen(sol.Ca);
	Cbeta = example_cast_eigen(sol.Cb);
	Calpha_occ = example_cast_eigen(Ca_occ);
	Cbeta_occ = example_cast_eigen(Cb_occ);

	 /// Orbital energies
	//  arma::vec Ea, Eb;

	  /// Fock operators
	Halpha = example_cast_eigen(sol.Ha);
	Hbeta = example_cast_eigen(sol.Hb);
	  /// Density matrices
	Dalpha = example_cast_eigen(sol.Pa);
	Dbeta = example_cast_eigen(sol.Pb);
	//P = sol.Pa + sol.Pb;

	  /// Coulomb operator
	J = example_cast_eigen(sol.J);
	  /// Exchange operators
	Kalpha = example_cast_eigen(sol.Ka);
	Kbeta = example_cast_eigen(sol.Kb);
	  /// KS-XC matrix
	XCalpha = example_cast_eigen(sol.XCa);
	XCbeta = example_cast_eigen(sol.XCb);

	  /// Imaginary exchange
	  //Matrix Ka_im, Kb_im;

	  /// Complex orbitals (for SIC)
	  //Matrix cCa, cCb;
	  /// Imaginary part of complex-CMO density matrix (for complex exchange contribution)
	  // Matrix Pa_im, Pb_im;

	  /// Energy information
	  //energy_t en;
	//sol.en = Ea + Eb;
	//Ea = sol.en/2;
	//Eb = sol.en/2;
	cout << "END wrap_force_UDFT\n";

	//return F2Erk;
	return G1Erk;
}







// FROM https://stackoverflow.com/questions/46700560/converting-an-armadillo-matrix-to-an-eigen-matrixd-and-vice-versa
// Converting an Armadillo Matrix to an Eigen MatriXd and vice versa
//One quick remark: If you are using Eigen's Mapping function, then you should automatically have the change in the Armadillo matrix (and vice versa), e.g.
//#include <RcppArmadillo.h>
//#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
//Eigen::MatrixXd HartreeFockClass::example_cast_eigen(arma::mat arma_A) {
Matrix HartreeFockClass::example_cast_eigen(arma::mat arma_A) {

  Eigen::MatrixXd eigen_B = Eigen::Map<Eigen::MatrixXd>(arma_A.memptr(),
                                                        arma_A.n_rows,
                                                        arma_A.n_cols);

  return eigen_B;
}

// [[Rcpp::export]]
//arma::mat HartreeFockClass::example_cast_arma(Eigen::MatrixXd eigen_A) {

arma::mat HartreeFockClass::example_cast_arma(Matrix eigen_A) {

  arma::mat arma_B = arma::mat(eigen_A.data(), eigen_A.rows(), eigen_A.cols(),
                               true, true); // ORIGINAL CODE IS  /*copy_aux_mem*/false, /*strict*/false);

  return arma_B;
}
/*
arma::mat HartreeFockClass::example_cast_arma(Matrix eigen_A) {


  //auto B[eigen_A.rows()][eigen_A.cols()] = eigen_A.data();
  auto B = eigen_A.data();
  //cout<<"\n Printing Matrix : \n";
  printDoubleMatrix(B,eigen_A.rows(),eigen_A.cols());

  arma::mat arma_B = arma::mat(B, eigen_A.rows(), eigen_A.cols(),
                               true, true); // ORIGINAL CODE IS  //copy_aux_mem//false, //strict//false);

  return arma_B;
}
*/
/*
arma::mat HartreeFockClass::example_cast_arma(Eigen::MatrixXcd eigen_A) {

  arma::mat arma_B = arma::mat(eigen_A.data(), eigen_A.rows(), eigen_A.cols(),
                               false, false);

  return arma_B;
}
*/
/*
void HartreeFockClass::printStdVector(std::vector<double> vec) {
	for (std::vector<int>::size_type i = 0; i < vec.size(); i++) {
    	cout <<' = (';
      	cout << vec.at(i) << ' ';
      	cout <<')\n';
    }
}
*/













//**R
//(x = matrix(1:4, ncol = 2))
//example_cast_eigen(x)
//example_cast_arma(x)

/*
//dft_t parse_dft(const Settings & set, bool init)
dft_t HartreeFockClass::parse_dft(bool init) {
  dft_t dft;
  dft.gridtol=0.0;
  dft.nl=false;
  dft.vv10_b=0.0;
  dft.vv10_C=0.0;

  // Use Lobatto quadrature?
  //dft.lobatto=set.get_bool("DFTLobatto");
  dft.lobatto=false;

  // Tolerance
  std::string tolkw = init ? "DFTInitialTol" : "DFTFinalTol";

  // Use static grid?
  //if(stricmp(set.get_string("DFTGrid"),"Auto")!=0) {
  //if(stricmp("50 -194","Auto")!=0) {
  if(stricmp("Auto","Auto")!=0) {
    //std::vector<std::string> opts=splitline(set.get_string("DFTGrid"));
	std::vector<std::string> opts=splitline("50 -194");
    if(opts.size()!=2) {
      throw std::runtime_error("Invalid DFT grid specified.\n");
    }

    dft.adaptive=false;
    dft.nrad=readint(opts[0]);
    dft.lmax=readint(opts[1]);
    if(dft.nrad<1 || dft.lmax==0) {
      throw std::runtime_error("Invalid DFT radial grid specified.\n");
    }

    // Check if l was given in number of points
    if(dft.lmax<0) {
      // Try to find corresponding Lebedev grid
      for(size_t i=0;i<sizeof(lebedev_degrees)/sizeof(lebedev_degrees[0]);i++)
	if(lebedev_degrees[i]==-dft.lmax) {
	  dft.lmax=lebedev_orders[i];
	  break;
	}
      if(dft.lmax<0) {
	std::ostringstream oss;
	oss << "Invalid DFT angular grid specified. Supported Lebedev grids:";
	for(size_t i=0;i<sizeof(lebedev_degrees)/sizeof(lebedev_degrees[0]);i++) {
	  if(i)
	    oss << ",";
	  oss << " " << lebedev_degrees[i];
	}
	oss << ".\n";
	throw std::runtime_error(oss.str());
      }
    }

  } else {
    dft.adaptive=true;
    dft.gridtol=set.get_double(tolkw);
  }

  // Parse functionals
  //parse_xc_func(dft.x_func,dft.c_func,set.get_string("Method"));
  //parse_xc_func(dft.x_func,dft.c_func,"UDFT");
  parse_xc_func(dft.x_func,dft.c_func,"LDA-LDA");

  // Parse VV10
  //std::string vv10s(set.get_string("VV10"));
  std::string vv10s("Auto");
  if(stricmp(vv10s,"Auto")==0) {
    // Determine automatically if VV10 is necessary
    if(dft.x_func>0)
      dft.nl=needs_VV10(dft.x_func,dft.vv10_b,dft.vv10_C);
    if(!dft.nl && dft.c_func>0)
      dft.nl=needs_VV10(dft.c_func,dft.vv10_b,dft.vv10_C);

  } else if(stricmp(vv10s,"True")==0 || stricmp(vv10s,"Yes")==0) {
    dft.nl=true;

    std::vector<std::string> vvopts=splitline(set.get_string("VV10Pars"));
    if(vvopts.size()!=2)
      throw std::runtime_error("Invalid VV10Pars!\n");

    dft.vv10_b=readdouble(vvopts[0]);
    dft.vv10_C=readdouble(vvopts[1]);

  } else if(stricmp(vv10s,"False")==0 || stricmp(vv10s,"No")==0) {
    // Do nothing

  } else if(vv10s.size()) {

    throw std::runtime_error("Error parsing VV10 setting.\n");
  }

  if(dft.nl) {
    if(dft.vv10_b <= 0.0 || dft.vv10_C <= 0.0) {
      std::ostringstream oss;
      oss << "VV10 parameters given b = " << dft.vv10_b << ", C = " << dft.vv10_C << " are not valid.\n";
      throw std::runtime_error(oss.str());
    }

    if(dft.adaptive)
      throw std::runtime_error("Adaptive DFT grids not supported with VV10.\n");

    std::vector<std::string> opts=splitline(set.get_string("NLGrid"));
    dft.nlnrad=readint(opts[0]);
    dft.nllmax=readint(opts[1]);
    if(dft.nlnrad<1 || dft.nllmax==0) {
      throw std::runtime_error("Invalid DFT radial grid specified.\n");
    }

    // Check if l was given in number of points
    if(dft.nllmax<0) {
      // Try to find corresponding Lebedev grid
      for(size_t i=0;i<sizeof(lebedev_degrees)/sizeof(lebedev_degrees[0]);i++)
	if(lebedev_degrees[i]==-dft.nllmax) {
	  dft.nllmax=lebedev_orders[i];
	  break;
	}
      if(dft.nllmax<0)
	throw std::runtime_error("Invalid DFT angular grid specified.\n");
    }

    // Check that xc grid is larger than nl grid
    if(dft.nrad < dft.nlnrad || dft.lmax < dft.nllmax)
      throw std::runtime_error("xc grid should be bigger than nl grid!\n");
  }

  return dft;
}
*/

















//Matrix HartreeFockClass::compute_2body_fock_general(const BasisSet& obs, const Matrix& D,
//        const BasisSet& D_bs, bool D_is_shelldiagonal,
//        double precision) {



//# 140 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 147 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/global.h"
//# 18 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in" 2



/*
//# 30 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
//void SCF::Fock_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid) const
//void SCF::Fock_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid) const
void HartreeFockClass::Fock_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid) const
//# 44 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
{
  //Timer t;
  //unsigned long int Nbf =	obs.nbf();

  unsigned long int Nbf = sol.Ha.n_cols;

  if(sol.P.n_rows != Nbf || sol.P.n_cols != Nbf) {
	  std::ostringstream oss;
  	  oss << "P" << " should be " << Nbf << " x " << Nbf << " but is " << sol.P.n_rows << " x " << sol.P.n_cols << "!\n";
  	  throw std::runtime_error(oss.str());
  };
  if(sol.Pa.n_rows != Nbf || sol.Pa.n_cols != Nbf) {
	  std::ostringstream oss;
	  oss << "Pa" << " should be " << Nbf << " x " << Nbf << " but is " << sol.Pa.n_rows << " x " << sol.Pa.n_cols << "!\n";
	  throw std::runtime_error(oss.str());
  };
  if(sol.Pb.n_rows != Nbf || sol.Pb.n_cols != Nbf) {
	  std::ostringstream oss;
	  oss << "Pb" << " should be " << Nbf << " x " << Nbf << " but is " << sol.Pb.n_rows << " x " << sol.Pb.n_cols << "!\n";
	  throw std::runtime_error(oss.str());
  };


  sol.J.zeros(Nbf,Nbf);

  sol.Ka.zeros(Nbf,Nbf);
  sol.Kb.zeros(Nbf,Nbf);




  sol.XCa.zeros(Nbf,Nbf);
  sol.XCb.zeros(Nbf,Nbf);


  double omega, kfull, kshort;
  range_separation(dft.x_func,omega,kfull,kshort);


//# 87 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  arma::mat Kafull, Kafull_im, Kbfull, Kbfull_im, Kashort, Kashort_im, Kbshort, Kbshort_im;
  //Matrix Kafull, Kbfull, Kashort, Kbshort;
  //Matrix Kafull_im, Kbfull_im, Kashort_im, Kbshort_im;
  Kafull.zeros(Nbf,Nbf);
  Kafull_im.zeros(Nbf,Nbf);
  Kbfull.zeros(Nbf,Nbf);
  Kbfull_im.zeros(Nbf,Nbf);
  Kashort.zeros(Nbf,Nbf);
  Kashort_im.zeros(Nbf,Nbf);
  Kbshort.zeros(Nbf,Nbf);
  Kbshort_im.zeros(Nbf,Nbf);


  memset(&en, 0, sizeof(energy_t));

  if(densityfit) {
    if(verbose) {
      printf("Forming density fitted Coulomb matrix ... ");
      fflush(stdout);
      t.set();
    }
    J=dfit.calcJ(sol.P);
    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
      fflush(stdout);
    }

    if(kfull!=0.0) {
      if(verbose) {
	printf("Forming density fitted exchange matrix ... ");
	fflush(stdout);
	t.set();
      }
//# 134 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
      if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {

	MatrixXcf cKa(dfit.calcK(sol.cCa,occa,fitmem));
	//Kafull=arma::real(cKa);
	Kafull=cKa.real();
	//Kafull_im=arma::imag(cKa);
	Kafull_im=cKa.imag();
      } else
	Kafull=dfit.calcK(sol.Ca,occa,fitmem);
      if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {

	MatrixXcf cKb(dfit.calcK(sol.cCb,occb,fitmem));
	Kbfull=cKb.real();
	Kbfull_im=cKb.imag());
      } else
	Kbfull=dfit.calcK(sol.Cb,occb,fitmem);

      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }

      if(kshort!=0.0) {
	if(verbose) {
	  printf("Forming density fitted short-range exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}
//# 173 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	  MatrixXcf cKa(dfit_rs.calcK(cCa,occa,fitmem));
	  Kashort=cKa.real();
	  Kashort_im=cKa.real();
	} else
	  Kashort=dfit_rs.calcK(sol.Ca,occa,fitmem);
	if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	  MatrixXcf cKb(dfit_rs.calcK(sol.cCb,occb,fitmem));
	  Kbshort=cKb.real();
	  Kbshort_im=cKb.imag();
	} else
	  Kbshort=dfit_rs.calcK(sol.Cb,occb,fitmem);

	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}
      }


    }

  } else {
    if(cholesky) {
      if(verbose) {
	printf("Forming Cholesky Coulomb matrix ... ");
	fflush(stdout);
	t.set();
      }
      J=chol.calcJ(sol.P);
      if(verbose) {
	printf("done (%s)\n",t.elapsed().c_str());
	fflush(stdout);
      }

      if(kfull!=0.0) {

	if(verbose) {
	  printf("Forming Cholesky exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}
//# 234 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	  MatrixXcf cKa(chol.calcK(sol.cCa,occa));
	  Kafull=cKa.real();
	  Kafull_im=cKa.imag();
	} else
	  Kafull=chol.calcK(sol.Ca,occa);
	if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	  MatrixXcf cKb(chol.calcK(sol.cCb,occb));
	  Kbfull=cKb.real();
	  Kbfull_im=cKb.imag();
	} else
	  Kbfull=chol.calcK(sol.Cb,occb);

	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}
      }

      if(kshort!=0.0) {
	if(verbose) {
	  printf("Forming short-range Cholesky exchange matrix ... ");
	  fflush(stdout);
	  t.set();
	}
//# 274 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {
	  MatrixXcf cKa(chol_rs.calcK(sol.cCa,occa));
	  Kashort=cKa.real();
	  Kashort_im=cKa.imag();
	} else
	  Kashort=chol_rs.calcK(sol.Ca,occa);
	if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {
	  MatrixXcf cKb(chol_rs.calcK(sol.cCb,occb));
	  Kbshort=cKb.real();
	  Kbshort_im=cKb.imag();
	} else
	  Kbshort=chol_rs.calcK(sol.Cb,occb);

	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}
      }

    } else {
      if(direct) {
//# 419 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

	if(kfull!=0.0) {
	  if(verbose) {
	    printf("Computing HF Coulomb and exchange matrices.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {

	    MatrixXcf cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    MatrixXcf cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    MatrixXcf cKa, cKb;
	    if(!decfock) {
	      scr.calcJK(cPa,cPb,J,cKa,cKb,intthr);
	    } else {

	      MatrixXcf Pahlp=decconv*cPa*decconv.trans();
	      MatrixXcf Pbhlp=decconv*cPb*decconv.trans();

	      scr.calcJK(Pahlp,Pbhlp,J,cKa,cKb,intthr);

	      J=decconv.trans()*J*decconv;
	      cKa=decconv.trans()*cKa*decconv;
	      cKb=decconv.trans()*cKb*decconv;
	    }
	    Kafull=cKa.real();
	    Kafull_im=cKa.imag();
	    Kbfull=cKb.real();
	    Kbfull_im=cKb.imag();
	  } else {
	    if(!decfock) {
	      scr.calcJK(sol.Pa,sol.Pb,J,Kafull,Kbfull,intthr);
	    } else {

	      Matrix Pahlp=decconv*sol.Pa*decconv.trans();
	      Matrix Pbhlp=decconv*sol.Pb*decconv.trans();

	      scr.calcJK(Pahlp,Pbhlp,sol.J,Kafull,Kbfull,intthr);

	      //sol.J=decconv.trans()*J*decconv;
	      sol.J=decconv.transpose()*J*decconv;
	      Kafull=decconv.transpose()*Kafull*decconv;
	      Kbfull=decconv.transpose()*Kbfull*decconv;
	    }
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }

	} else {
	  if(verbose) {
	    printf("Computing HF Coulomb matrix.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(!decfock)
		  sol.J=scr.calcJ(sol.P,intthr);
	  else {
	    Matrix Phlp=decconv*sol.P*decconv.trans();
	    Matrix Jhlp=scr.calcJ(Phlp,intthr);
	    sol.J=decconv.transpose()*Jhlp*decconv;
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}

	if(kshort!=0.0) {
	  if(verbose) {
	    printf("Computing screened exchange matrix.\nScreening integrals with tolerance %.3e ... ",intthr);
	    fflush(stdout);
	    t.set();
	  }
	  if(Pa_im.n_rows == Pa.n_rows && Pa_im.n_cols == Pa.n_cols) {

	    MatrixXcf cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    MatrixXcf cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    MatrixXcf cKa, cKb;
	    if(!decfock) {
	      scr_rs.calcK(cPa,cPb,cKa,cKb,intthr);
	    } else {

	      MatrixXcf Pahlp=decconv*cPa*decconv.trans();
	      MatrixXcf Pbhlp=decconv*cPb*decconv.trans();

	      scr_rs.calcK(Pahlp,Pbhlp,cKa,cKb,intthr);

	      cKa=decconv.transpose()*cKa*decconv;
	      cKb=decconv.transpose()*cKb*decconv;
	    }
	    Kashort=cKa.real();
	    Kashort_im=cKa.imag();
	    Kbshort=cKb.real();
	    Kbshort_im=cKb.imag();
	  } else {
	    if(!decfock) {
	      scr_rs.calcK(sol.Pa,sol.Pb,Kashort,Kbshort,intthr);
	    } else {

	      Matrix Pahlp=decconv*sol.Pa*decconv.transpose();
	      Matrix Pbhlp=decconv*sol.Pb*decconv.transpose();

	      scr_rs.calcK(Pahlp,Pbhlp,Kashort,Kbshort,intthr);

	      Kashort=decconv.transpose()*Kashort*decconv;
	      Kbshort=decconv.transpose()*Kbshort*decconv;
	    }
	  }
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}

      } else {

	if(verbose) {
	  printf("Forming HF Coulomb matrix ... ");
	  fflush(stdout);
	  t.set();
	}
	sol.J=tab.calcJ(sol.P);
	if(verbose) {
	  printf("done (%s)\n",t.elapsed().c_str());
	  fflush(stdout);
	}

	if(kfull!=0.0) {

	  if(verbose) {
	    printf("Forming HF exchange matrix ... ");
	    fflush(stdout);
	    t.set();
	  }
//# 592 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {

	    MatrixXcf cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    MatrixXcf cKa(tab.calcK(cPa));
	    Kafull=cKa.real();
	    Kafull_im=cKa.imag();
	  } else
	    Kafull=tab.calcK(Pa);
	  if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {

	    MatrixXcf cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    MatrixXcf cKb(tab.calcK(cPb));
	    Kbfull=cKb.real();
	    Kbfull_im=cKb.imag();
	  } else
	    Kbfull=tab.calcK(Pb);

	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}

	if(kshort!=0.0) {
	  if(verbose) {
	    printf("Forming short-range HF exchange matrix ... ");
	    fflush(stdout);
	    t.set();
	  }
//# 638 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
	  if(sol.Pa_im.n_rows == sol.Pa.n_rows && sol.Pa_im.n_cols == sol.Pa.n_cols) {

	    MatrixXcf cPa(sol.Pa*std::complex<double>(1.0,0.0) + sol.Pa_im*std::complex<double>(0.0,1.0));
	    MatrixXcf cKa(tab_rs.calcK(cPa));
	    Kashort=cKa.real();
	    Kashort_im=cKa.imag();
	  } else
	    Kashort=tab_rs.calcK(sol.Pa);
	  if(sol.Pb_im.n_rows == sol.Pb.n_rows && sol.Pb_im.n_cols == sol.Pb.n_cols) {

	    MatrixXcf cPb(sol.Pb*std::complex<double>(1.0,0.0) + sol.Pb_im*std::complex<double>(0.0,1.0));
	    MatrixXcf cKb(tab_rs.calcK(cPb));
	    Kbshort=cKb.real();
	    Kbshort_im=cKb.imag();
	  } else
	    Kbshort=tab_rs.calcK(sol.Pb);
	  if(verbose) {
	    printf("done (%s)\n",t.elapsed().c_str());
	    fflush(stdout);
	  }
	}
      }
    }
  }


//# 676 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  sol.Ka=kfull*Kafull + kshort*Kashort;
  if(sol.Pa.n_rows == sol.Pa_im.n_rows && sol.Pa.n_cols == sol.Pa_im.n_cols)
    Ka_im=kfull*Kafull_im + kshort*Kashort_im;
  else
	  sol.Ka_im.clear();
  sol.Kb=kfull*Kbfull + kshort*Kbshort;
  if(sol.Pb.n_rows == sol.Pb_im.n_rows && sol.Pb.n_cols == sol.Pb_im.n_cols)
	  sol.Kb_im=kfull*Kbfull_im + kshort*Kbshort_im;
  else
	  sol.Kb_im.clear();
//# 697 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

  en.Exc=0.0;
  if(dft.x_func>0 || dft.c_func>0) {
    if(verbose) {
      printf("Computing DFT exchange and correlation ... ");
      fflush(stdout);
      t.set();
    }
    double Nelnum;

    grid.eval_Fxc(dft.x_func,dft.c_func,sol.Pa,sol.Pb,sol.XCa,sol.XCb,en.Exc,Nelnum);



    double rel_diff=(Nelnum-Nel)*100.0/Nel;

    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
      printf("Numerically integrated density is %.5f (%+.4f %%).\n",Nelnum,rel_diff);
    }
    if(fabs(rel_diff)>1e-2) {
      std::ostringstream oss;

      oss << "Warning - numerically integrated density seems inaccurate.\n";
      if(verbose)
	cout << oss.str();

    }
  }

  en.Esic=0.0;

  if(dft.nl) {
    if(verbose) {
      printf("Computing non-local correlation ... ");
      fflush(stdout);
      t.set();
    }


    Matrix XC(sol.XCa);
    XC.zeros();
    grid.eval_VV10(nlgrid,dft.vv10_b,dft.vv10_C,P,XC,en.Enl);
    sol.XCa+=XC;
    sol.XCb+=XC;

    if(verbose) {
      printf("done (%s)\n",t.elapsed().c_str());
    }
  }


//# 785 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

  sol.Ha=Hcore+sol.J+sol.XCa;
  sol.Hb=Hcore+sol.J+sol.XCb;
  if(kshort!=0.0 || kfull!=0.0) {
	  sol.Ha-=sol.Ka;
	  sol.Hb-=sol.Kb;
    //en.Exc-=0.5*(arma::trace(sol.Pa*sol.Ka)+arma::trace(sol.Pb*sol.Kb));
	  en.Exc-=0.5*((sol.Pa.cwiseProduct(sol.Ka))+(sol.Pb.cwiseProduct(sol.Kb)));
    if(sol.Pa.n_rows == sol.Pa_im.n_rows && sol.Pa.n_cols == sol.Pa_im.n_cols)
      //en.Exc+=0.5*(arma::trace(sol.Pa_im*sol.Ka_im)+arma::trace(sol.Pb_im*sol.Kb_im));
    	en.Exc+=0.5*((sol.Pa_im.cwiseProduct(sol.Ka_im))+(sol.Pb_im.cwiseProduct(sol.Kb_im));
  }
//# 811 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

//# 822 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"

  if(freeze.size()>0) {

    freeze_orbs(freeze,sol.Ca,S,sol.Ha,verbose);
    freeze_orbs(freeze,sol.Cb,S,sol.Hb,verbose);
  }
  fflush(stdout);
  if(dimcalc) {

    std::vector<nucleus_t> nuclei(basisp->get_nuclei());
    for(size_t i=0;i<nuclei.size();i++)
      if(nuclei[i].r.x!=0.0 || nuclei[i].r.y!=0.0)
        throw std::logic_error("Nuclei must be on z axis for dimer calculation!\n");

    //arma::ivec mvals(basisp->get_m_values());
    Eigen::VectorXf mvals(basisp->get_m_values());
    mvals.sort();


    sol.Ha=block_m(sol.Ha,mvals);
    sol.Hb=block_m(sol.Hb,mvals);
  }

  //en.Ekin=arma::trace(sol.P*T);
  en.Ekin=sol.P.cwiseProduct(T);
  //en.Enuca=arma::trace(sol.P*Vnuc);
  en.Enuca=sol.P.cwiseProduct(Vnuc);
  en.Enucr=Enuc;
  //en.Eone=arma::trace(sol.P*Hcore);
  en.Eone=sol.P.cwiseProduct(Hcore);
  //en.Ecoul=0.5*arma::trace(sol.P*sol.J);
  en.Ecoul=sol.P.cwiseProduct(sol.J);

  en.Eel=en.Ecoul+en.Exc+en.Eone+en.Enl;
  en.E=en.Eel+en.Enucr;

//# 876 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
  //if(!arma::is_finite(sol.Ha)) {
  if(!Eigen::isfinite(sol.Ha)) {
//# 887 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Alpha Fock operator is not finite.\n");
  }
  if(!Eigen::isfinite(sol.Hb)) {
//# 900 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    throw std::runtime_error("Beta Fock operator is not finite.\n");
  }

  if(!std::isfinite(en.E)) {
//# 935 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-fock.cpp.in"
    std::ostringstream oss;
    oss << "\nSomething wrong with total energy " << en.E <<"?\nEnding program.\n";
    throw std::runtime_error(oss.str());
  }
}
*/



//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
//# 1 "/usr/include/stdc-predef.h" 1






//# 43 "/usr/include/stdc-predef.h"
//# 51 "/usr/include/stdc-predef.h"


//# 1 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in" 2



/*


//# 25 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
//arma::vec SCF::force_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol)
//std::vector<float> HartreeFockClass::force_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol)
arma::vec HartreeFockClass::force_UDFT(SCF ErkaleSCF, uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol)
//# 33 "/home/carlos/eclipse-workspace/2018-11-13-erkale/erkale-master/build/erkale-eclipse/src/scf-force.cpp.in"
{
  //arma::mat W;
  Matrix W;
  W=form_density(sol.Ea,sol.Ca,occa)+form_density(sol.Eb,sol.Cb,occb);

  //arma::vec fpul_kin=basisp->kinetic_pulay(sol.P);
  Eigen::VectorXf fpul_kin=basisp->kinetic_pulay(sol.P);

  Eigen::VectorXf fpul_nuc=basisp->nuclear_pulay(sol.P);

  Eigen::VectorXf fnuc=basisp->nuclear_der(sol.P);

  Eigen::VectorXf forth=basisp->overlap_der(W);

  Eigen::VectorXf frep=basisp->nuclear_force();


  double omega, kfull, kshort;
  range_separation(dft.x_func,omega,kfull,kshort);

  Eigen::VectorXf fx_full;
  fx_full.zeros(fnuc.n_elem);

  Eigen::VectorXf fx_short;
  fx_short.Zero(fnuc.n_elem);
  if(kfull==0.0) {
    if(densityfit) {
      if(kshort!=0.0 || kfull!=0.0)
	throw std::runtime_error("Forces not implemented for density fitting of exact exchange.\n");
      fx_full=dfit.forceJ(sol.P);
    } else {
      if(!direct)
	scr.fill(basisp,intthr,verbose);
      fx_full=scr.forceJ(sol.P,tol);
    }
  } else {
    if(!direct)
      scr.fill(basisp,intthr,verbose);
    fx_full=scr.forceJK(sol.Pa,sol.Pb,tol,kfull);
  }
  if(omega != 0.0) {

    scr_rs.set_range_separation(omega,0.0,1.0);
    scr_rs.fill(basisp,intthr,verbose);

    fx_short=scr_rs.forceK(sol.Pa,sol.Pb,tol,kshort);
  }

  //arma::vec fxc;
  Eigen::VectorXf fxc;
  fxc.Zero(fnuc.n_elem);
  if(dft.x_func>0 || dft.c_func>0) {
    fxc=grid.eval_force(dft.x_func,dft.c_func,sol.Pa,sol.Pb);
  }

  if(dft.nl) {
	Eigen::VectorXf vv10f(grid.eval_VV10_force(nlgrid,dft.vv10_b,dft.vv10_C,sol.P));

    fxc+=vv10f;
  }
  Eigen::VectorXf ftot=fpul_kin+fpul_nuc+fnuc+forth+frep;
  Eigen::VectorXf ffull=ftot+fx_full;

  ffull+=fxc+fx_short;



  return ffull;
}

*/




















































































































//****************************************************************************************************************************************
//************************************ END DFT PART FROM ERKALE **************************************************************************
//****************************************************************************************************************************************


} /* namespace HartreeFock */

using namespace HartreeFock;
//using HartreeFock::HartreeFockClass;


// REMOVE MAIN TO BE ALLOWED TO INCLUDE THIS CLASS IN THE LIBRARY. THE EXECUTABLE WIL BE BUILT WITH main.cpp INSTEAD OF THIS.
/*
int main(int argc, char* argv[]) {
	HartreeFockClass hf;
	hf.run("h2o.xyz","aug-cc-pVDZ");



	cout <<"hello world\n";

		size_t N = 3;
		vector<double> A(N * N), B(N * N), C(N * N);

//		HartreeFock::HartreeFockClass hf;

		//Fill A and B with random numbers in the range [0,1]:
		hf.rnd_fill(A, 0.0, 1.0);
		hf.rnd_fill(B, 0.0, 1.0);

		//Multiply matrices A and B, save the result in C
		int sz = N;
		//hf.matrix_multiply(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
		//hf.mat_mult(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
		//hf.mat_mult(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
		//hf.mat_mult(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
		//hf.mat_mult(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
		hf.mat_mult(A, sz, sz, B, sz, sz, C, sz, sz);


		//print A, B, C on the standard output
		hf.print_matrix(A, sz, sz);
		hf.print_matrix(B, sz, sz);
		hf.print_matrix(C, sz, sz);


		return 0;


	return 0;
}
*/
