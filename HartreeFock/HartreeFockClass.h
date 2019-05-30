/*
 * HartreeFockClass.h
 *
 *  Created on: 18 mar. 2017
 *      Author: carlos
 */

// ADDED TO AVOID DOUBLE INCLUDE THIS HEADER FILE.
//#ifndef INCLUDED_HARTREEFOCKCLASS_H
//#define INCLUDED_HARTREEFOCKCLASS_H
// END ADDED TO AVOID DOUBLE INCLUDE THIS HEADER FILE.

#ifndef HARTREEFOCKCLASS_H_
#define HARTREEFOCKCLASS_H_

// standard C++ headers
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

// Eigen matrix algebra library
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
//#include <Eigen/IOFormat> // IS A CLASS.

// have BTAS library?
#ifdef LIBINT2_HAVE_BTAS
#include <btas/btas.h>
#endif  // LIBINT2_HAVE_BTAS

// Libint Gaussian integrals library
#include <libint2/diis.h>
#include <libint2/util/intpart_iter.h>
#include <libint2/chemistry/sto3g_atomic_density.h>
#include <libint2.hpp>

#if defined(_OPENMP)
#include <omp.h>
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
























// uncomment if want to report integral timings
// N.B. integral engine timings are controled in engine.h
#define REPORT_INTEGRAL_TIMINGS

//typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix;  // import dense, dynamically sized Matrix type from Eigen;
             // this is a matrix with row-major storage
             // (http://en.wikipedia.org/wiki/Row-major_order)
// to meet the layout of the integrals returned by the Libint integral library
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic>
    DiagonalMatrix;


// TO USE WITH DFT FROM ERKALE HELGAKER COMPLEX MATRICES.
//typedef Eigen::MatrixXcf<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
typedef Eigen::Matrix<std::complex< double >, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    MatrixXcf;  // import dense, dynamically sized Matrix type from Eigen;
             // this is a matrix with row-major storage
             // (http://en.wikipedia.org/wiki/Row-major_order)





//using DFTGrid;


using namespace std;
using std::vector;

using libint2::Shell;
using libint2::Atom;
using libint2::BasisSet;
using libint2::Operator;
using libint2::BraKet;

using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;

typedef std::vector<Matrix> VectorOfMatrices;
/*
// OPERATOR IMPLEMENTATION TO PERFORM SUMMATION OF 2 STD VECTORS AS DONE BY +.
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());// OPERATOR IMPLEMENTATION TO PERFORM SUMMATION OF 2 STD VECTORS AS DONE BY +.

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<T>());
    return result;
}
*/


namespace libint2 {
// ADDED TO AVOID DOUBLE DEFINITION OF libint2::nthreads IN main.o.
#ifndef DEFINED_LIBINT2NTHREADS_
#define DEFINED_LIBINT2NTHREADS_
// int nthreads; // ORIGINAL CODE, DEFINITION, CAUSES DOUBLE DEFINITION WHEN USING main.cpp.
extern int nthreads; // DECLARATION INSTEAD OF DEFINITION.
#endif

/// fires off \c nthreads instances of lambda in parallel
template <typename Lambda>
void parallel_do(Lambda& lambda) {
#ifdef _OPENMP
#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#else  // use C++11 threads
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id != libint2::nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(lambda, thread_id));
    else
      lambda(thread_id);
  }  // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
#endif
}
}



namespace HartreeFock {


//Fortran subroutine definition "translated" to C++
/* 2018-04-30 FROM FORTRAN MODULE MATRIXMULTIPLY, ELIMINATED.
extern "C" {
    void matrix_multiply(double *A, int *rows_A, int *cols_A, double *B, int *rows_B, int *cols_B, double *C, int *rows_C, int *cols_C);

//    void ANTC(bool *UHF, unsigned long *JCycle, unsigned long *inputNCycles,char *inputjobname, unsigned long *inputjobname_len, double *D,double *pivHam, double *pivFock, double *pivCoulomb, double *pivExchange, double *pivOver,
//    		unsigned long *NBasis, unsigned long *inputNSpin, unsigned long *inputNAtoms, unsigned long *inputNShell,long *inputJAN, unsigned long *inputAOS, unsigned long *inputShellT, unsigned long *inputShellC,
//			unsigned long *inputNAE, unsigned long *inputNBE, unsigned long *inputNE, unsigned long *inputIAN, double *inputAtmCo,double *inputAtmChg,
//            double *denerrj, double *Crit,bool *ANTOn);

};
*/

//Fortran subroutine definition "translated" to C++

extern "C" {
//    void ANTC(bool *UHF, int *JCycle, int *inputNCycles,char *inputjobname, int *inputjobname_len,double *D,double *pivHam,double *pivFock,double *pivCoulomb,double *pivExchange,double *pivOver,
//    		int *NBasis,int *inputNSpin,int *inputNAtoms,int *inputNShell,int *inputJAN,int *inputAOS,int *inputShellT,int *inputShellC,
//			double *inputNAE,double *inputNBE,double *inputNE,int *inputIAN,double *inputAtmCo,double *inputAtmChg,
 //           double *denerrj,double *Crit,bool *ANTOn);

//      void antc(bool *UHF, unsigned long *JCycle, unsigned long *inputNCycles,char *inputjobname, unsigned long *inputjobname_len, double *D,double *pivHam, double *pivFock, double *pivCoulomb, double *pivExchange, double *pivOver,
//    		unsigned long *NBasis, unsigned long *inputNSpin, unsigned long *inputNAtoms, unsigned long *inputNShell,long *inputJAN, unsigned long *inputAOS, unsigned long *inputShellT, unsigned long *inputShellC,
//			unsigned long *inputNAE, unsigned long *inputNBE, unsigned long *inputNE, unsigned long *inputIAN, double *inputAtmCo,double *inputAtmChg,
//            double *denerrj, double *Crit,bool *ANTOn);
        void antc(bool *UHF,unsigned long *JCycle,unsigned long *inputNCycles,char *inputjobname,unsigned long *inputjobname_len,//);
        		  unsigned long *NBasis,unsigned long *inputNSpin,unsigned long *inputNAtoms,unsigned long *inputNShell,
                  double *DA,double *DB,
				  double *pivHam,
				  double *pivFockA,double *pivFockB,
				  double *pivCoulomb,
				  double *pivExchangeA,double *pivExchangeB,
				  double *pivOver,
				  double *outFockHWA,double *outFockHWB,
				  double *outGibbsYA,double *outGibbsYB,double *outGibbsYKernel1A,double *outGibbsYKernel1B,double *outGibbsYKernel2A,double *outGibbsYKernel2B,
				  unsigned long *inputIAN,long *inputJAN,unsigned long *inputAOS,unsigned long *inputShellT,unsigned long *inputShellC,
				  unsigned long *inputNAE,unsigned long *inputNBE,unsigned long *inputNE,double *inputAtmCo,double *inputAtmChg,
				  double *denerrj,double *Crit,bool *ANTOn);
};


//extern "C" typedef void matrix_multiply(double *A, int *rows_A, int *cols_A, double *B, int *rows_B, int *cols_B, double *C, int *rows_C, int *cols_C);

class HartreeFockClass {
public:
	HartreeFockClass();
	virtual ~HartreeFockClass();

	void rnd_fill(vector<double> &V, double lower, double upper);
	void print_matrix(vector<double> const &V, int rows, int cols);

	void writeEigenMatrixToCSVfile(string name, Matrix matrix);
	//void writeEigenMatrixToCSVfile(const char * name, Matrix matrix);

	//friend void matrix_multiply(double *A, int *rows_A, int *cols_A, double *B, int *rows_B, int *cols_B, double *C, int *rows_C, int *cols_C);
	//matrix_multiply mat_mult;
	//void mat_mult(std::vector<double, std::allocator<double>> A, int rows_A, int cols_A, std::vector<double, std::allocator<double>> B, int rows_B, int cols_B, std::vector<double, std::allocator<double>> C, int rows_C, int cols_C);
	//void mat_mult(vector<double> const &A, int rows_A, int cols_A,vector<double> const &B, int rows_B, int cols_B,vector<double> const &C, int rows_C, int cols_C);

	/* 2018-04-30 FROM FORTRAN MODULE MATRIXMULTIPLY, ELIMINATED.
	void mat_mult(vector<double> &A, int rows_A, int cols_A,vector<double> &B, int rows_B, int cols_B,vector<double> &C, int rows_C, int cols_C);
    */

	void run(char *geom, char *basis);

	void runant(const char *geomname, const char *basisname, const char *fullname, const bool debugflag, const bool erkaleflag);


	int **alloc_2d_int(int rows, int cols);
	double **alloc_2d_double(int rows, int cols);


	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        	Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 	// this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 	// to meet the layout of the integrals returned by the Libint integral library
	typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic>
        	DiagonalMatrix;






//-------------------------------------------------------------------------------------------------------------
//----------------- DEFINITIONS FROM HARTREE-FOCK -------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//void printDoubleMatrix(double ** mat, int ROW, int COL);
void printDoubleMatrix(double * mat, int ROW, int COL);
//void printDoubleMatrix(double mat[], int ROW, int COL);
void printEigenmatrix(Matrix &m1); // ADDED BY C.SALGADO.
void printEigenmatrixCommaInitFmt(Matrix &m1); // ADDED BY C.SALGADO.
void printEigenmatrixCleanFmt(Matrix &m1); // ADDED BY C.SALGADO.
void printEigenmatrixOctaveFmt(Matrix &m1); // ADDED BY C.SALGADO.
void printEigenmatrixHeavyFmt(Matrix &m1); // ADDED BY C.SALGADO.

std::vector<Atom> read_geometry(const std::string& filename);
Matrix compute_soad(const std::vector<Atom>& atoms);
// computes norm of shell-blocks of A
Matrix compute_shellblock_norm(const BasisSet& obs, const Matrix& A);

template <Operator obtype>
std::array<Matrix, libint2::operator_traits<obtype>::nopers> compute_1body_ints(
    const BasisSet& obs, const std::vector<Atom>& atoms = std::vector<Atom>());

#if LIBINT2_DERIV_ONEBODY_ORDER
template <Operator obtype>
std::vector<Matrix> compute_1body_ints_deriv(unsigned deriv_order,
                                             const BasisSet& obs,
                                             const std::vector<Atom>& atoms);
#endif  // LIBINT2_DERIV_ONEBODY_ORDERR

template <libint2::Operator Kernel = libint2::Operator::coulomb>
Matrix compute_schwarz_ints(
    const BasisSet& bs1, const BasisSet& bs2 = BasisSet(),
    bool use_2norm = false,  // use infty norm by default
    typename libint2::operator_traits<Kernel>::oper_params_type params =
        libint2::operator_traits<Kernel>::default_params());
Matrix compute_do_ints(const BasisSet& bs1, const BasisSet& bs2 = BasisSet(),
                       bool use_2norm = false  // use infty norm by default
                       );

//using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
shellpair_list_t obs_shellpair_list;  // shellpair list for OBS

/// computes non-negligible shell pair list; shells \c i and \c j form a
/// non-negligible
/// pair if they share a center or the Frobenius norm of their overlap is
/// greater than threshold
shellpair_list_t compute_shellpair_list(const BasisSet& bs1,
                                        const BasisSet& bs2 = BasisSet(),
                                        double threshold = 1e-12);

Matrix compute_2body_fock(
    const BasisSet& obs, const Matrix& D,
    double precision = std::numeric_limits<
        double>::epsilon(),  // discard contributions smaller than this
    const Matrix& Schwarz = Matrix()  // K_ij = sqrt(||(ij|ij)||_\infty); if
                                       // empty, do not Schwarz screen
    );
VectorOfMatrices compute_2body_fock_uhf(
    const BasisSet& obs, const Matrix& Dalpha, const Matrix& Dbeta,
    double precision = std::numeric_limits<
        double>::epsilon(),  // discard contributions smaller than this
    const Matrix& Schwarz = Matrix()  // K_ij = sqrt(||(ij|ij)||_\infty); if
                                       // empty, do not Schwarz screen
    );
VectorOfMatrices compute_2body_JK_uhf(
    const BasisSet& obs, const Matrix& Dalpha, const Matrix& Dbeta,
    double precision = std::numeric_limits<
        double>::epsilon(),  // discard contributions smaller than this
    const Matrix& Schwarz = Matrix()  // K_ij = sqrt(||(ij|ij)||_\infty); if
                                       // empty, do not Schwarz screen
    );
// an Fock builder that can accept densities expressed a separate basis
Matrix compute_2body_fock_general(
    const BasisSet& obs, const Matrix& D, const BasisSet& D_bs,
    bool D_is_sheldiagonal = false,  // set D_is_shelldiagonal if doing SOAD
    double precision = std::numeric_limits<
        double>::epsilon()  // discard contributions smaller than this
    );

#if LIBINT2_DERIV_ERI_ORDER
template <unsigned deriv_order>
std::vector<Matrix> compute_2body_fock_deriv(
    const BasisSet& obs, const std::vector<Atom>& atoms,
    const Matrix& D,
    double precision = std::numeric_limits<
        double>::epsilon(),  // discard contributions smaller than this
    const Matrix& Schwarz = Matrix()  // K_ij = sqrt(||(ij|ij)||_\infty); if
                                       // empty, do not Schwarz screen
    );
#endif  // LIBINT2_DERIV_ERI_ORDER

#if LIBINT2_DERIV_ERI_ORDER
template <unsigned deriv_order>
std::vector<VectorOfMatrices> compute_2body_fock_uhf_deriv(
//std::vector<std::vector<Matrix>> compute_2body_fock_uhf_deriv(
    const BasisSet& obs, const std::vector<Atom>& atoms,
    const Matrix& Dalpha, const Matrix& Dbeta,
    double precision = std::numeric_limits<
        double>::epsilon(),  // discard contributions smaller than this
    const Matrix& Schwarz = Matrix()  // K_ij = sqrt(||(ij|ij)||_\infty); if
                                       // empty, do not Schwarz screen
    );
#endif  // LIBINT2_DERIV_ERI_ORDER



/*
#if LIBINT2_DERIV_ERI_ORDER
template <unsigned deriv_order>
std::vector<VectorOfMatrices> compute_2body_JK_uhf_deriv(
//std::vector<std::vector<Matrix>> compute_2body_fock_uhf_deriv(
    const BasisSet& obs, const std::vector<Atom>& atoms,
    const Matrix& Dalpha, const Matrix& Dbeta,
    double precision = std::numeric_limits<
        double>::epsilon(),  // discard contributions smaller than this
    const Matrix& Schwarz = Matrix()  // K_ij = sqrt(||(ij|ij)||_\infty); if
                                       // empty, do not Schwarz screen
    );
#endif  // LIBINT2_DERIV_ERI_ORDER
*/




// returns {X,X^{-1},S_condition_number_after_conditioning}, where
// X is the generalized square-root-inverse such that X.transpose() * S * X = I
// columns of Xinv is the basis conditioned such that
// the condition number of its metric (Xinv.transpose . Xinv) <
// S_condition_number_threshold
std::tuple<Matrix, Matrix, double> conditioning_orthogonalizer(
    const Matrix& S, double S_condition_number_threshold);

#ifdef LIBINT2_HAVE_BTAS
#define HAVE_DENSITY_FITTING 1
struct DFFockEngine {
  const BasisSet& obs;
  const BasisSet& dfbs;
  DFFockEngine(const BasisSet& _obs, const BasisSet& _dfbs)
      : obs(_obs), dfbs(_dfbs) {}

  typedef btas::RangeNd<CblasRowMajor, std::array<long, 3>> Range3d;
  typedef btas::Tensor<double, Range3d> Tensor3d;
  Tensor3d xyK;

  // a DF-based builder, using coefficients of occupied MOs
  Matrix compute_2body_fock_dfC(const Matrix& Cocc);
};
#endif  // HAVE_DENSITY_FITTING


void api_basic_compile_test(const BasisSet& obs);






std::tuple<Matrix, Matrix, size_t, double, double> gensqrtinv(
    const Matrix& S, bool symmetric = false,
    double max_condition_number = 1e8);

//std::tuple<Matrix, Matrix, double> conditioning_orthogonalizer(
//    const Matrix& S, double S_condition_number_threshold);

Matrix compute_2body_2index_ints(const BasisSet& bs);





//-------------------------------------------------------------------------------------------------------------
//------------- END DEFINITIONS FROM HARTREE-FOCK -------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------

	//typedef libint2::Shell Shell;
	//typedef libint2::Atom Atom;
	//typedef libint2::BasisSet HBasisSet;
	//using libint2::Shell; // You can't bring a name from the class scope into namespace scope with using directive. You can only be using a name defined at namespace scope, in a different namespace.
	//using libint2::Atom; // You can't bring a name from the class scope into namespace scope with using directive. You can only be using a name defined at namespace scope, in a different namespace.
	//using libint2::BasisSet; // You can't bring a name from the class scope into namespace scope with using directive. You can only be using a name defined at namespace scope, in a different namespace.






















// DFT HELGAKER ERKALE USEFUL THINGS

/*
typedef struct {
  /// Orbitals
//  arma::mat Ca, Cb;
  Matrix Ca, Cb;
  /// Orbital energies
//  arma::vec Ea, Eb;
  std::vector<float> Ea, Eb;
  /// Fock operators
  Matrix Ha, Hb;
  /// Density matrices
  Matrix P, Pa, Pb;

  /// Coulomb operator
  Matrix J;
  /// Exchange operators
  Matrix Ka, Kb;
  /// KS-XC matrix
  Matrix XCa, XCb;

  /// Imaginary exchange
  Matrix Ka_im, Kb_im;

  /// Complex orbitals (for SIC)
  Matrix cCa, cCb;
  /// Imaginary part of complex-CMO density matrix (for complex exchange contribution)
  Matrix Pa_im, Pb_im;

  /// Energy information
  //energy_t en;
} uscf_t;
*/
/*
/// DFT settings
typedef struct {
  /// Used exchange functional
  int x_func;
  /// Used correlation functional
  int c_func;

  /// Adaptive grid?
  bool adaptive;
  /// Integration grid tolerance
  double gridtol;
  /// Amount of radial shells (if not adaptive)
  int nrad;
  /// Maximum angular quantum number to integrate exactly (if not adaptive)
  int lmax;
  /// Use Lobatto quadrature?
  bool lobatto;

  // Non-local part?
  bool nl;
  double vv10_b, vv10_C;
  int nlnrad, nllmax;
} dft_t;
*/


//  void Fock_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid) const;
//void Fock_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid);// const; COMMENTED BY C.SALGADO ON 2019-01-01.

//Eigen::VectorXf force_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol);
//  arma::vec force_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol);



SCF initialize_Erkale_SCF(string basename, ErkaleBasisSet& basisp, Settings& setp, unsigned long int Nbf, int Nel, int mult, Matrix& S, Matrix& T, Matrix& Vnuc, Matrix& Hcore );

//void wrap_Fock_UDFT(string basename, const double Ea, const double Eb, const int naocc, const int nbocc, const Matrix& Calpha, const Matrix& Cbeta,const Matrix& Calpha_occ, const Matrix& Cbeta_occ, const Matrix& Halpha, const Matrix& Hbeta, const Matrix& Dalpha, const Matrix& Dbeta, const Matrix& J, const Matrix& Kalpha, const Matrix& Kbeta, const Matrix& XCalpha, const Matrix& XCbeta);
void wrap_Fock_UDFT(SCF ErkaleSCF, ErkaleBasisSet& basis, Settings& set, string basename, unsigned long int Nbf, std::vector<double> Ea, std::vector<double> Eb, int naocc, int nbocc, Matrix& S, Matrix& Calpha, Matrix& Cbeta,Matrix& Calpha_occ, Matrix& Cbeta_occ, Matrix& Halpha, Matrix& Hbeta, Matrix& Dalpha, Matrix& Dbeta, Matrix& J, Matrix& Kalpha, Matrix& Kbeta, Matrix& XCalpha, Matrix& XCbeta);

std::vector<std::vector<Matrix> > wrap_force_UDFT(SCF ErkaleSCF, ErkaleBasisSet& basis, Settings& set, string basename, unsigned long int Nbf, std::vector<double> Ea, std::vector<double> Eb, int naocc, int nbocc, Matrix& S, Matrix& Calpha, Matrix& Cbeta, Matrix& Calpha_occ, Matrix& Cbeta_occ, Matrix& Halpha, Matrix& Hbeta, Matrix& Dalpha, Matrix& Dbeta, Matrix& J, Matrix& Kalpha, Matrix& Kbeta, Matrix& XCalpha, Matrix& XCbeta);

//Eigen::MatrixXd example_cast_eigen(arma::mat arma_A);
Matrix example_cast_eigen(arma::mat arma_A);

//arma::mat example_cast_arma(Eigen::MatrixXd eigen_A);
arma::mat example_cast_arma(Matrix eigen_A);
//arma::mat example_cast_arma(Eigen::MatrixXd eigen_A);
//arma::mat example_cast_arma(Eigen::MatrixXcd eigen_A);




template<typename T> inline void printStdVector(std::vector<T> vec){
	for (std::vector<int>::size_type i = 0; i < vec.size(); i++) {
    	cout <<" = (";
      	cout << vec.at(i) << ' ';
      	cout <<")\n";
    }
}
//dft_t parse_dft(bool init);



//void Fock_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid) const;

//Eigen::VectorXf force_UDFT(uscf_t & sol, const std::vector<double> & occa, const std::vector<double> & occb, const dft_t dft, DFTGrid & grid, DFTGrid & nlgrid, double tol);















































}; // END HartreeFockClass

#ifdef LIBINT2_HAVE_BTAS
# define HAVE_DENSITY_FITTING 1
	struct DFFockEngine {
		const BasisSet& obs;
		const BasisSet& dfbs;
		DFFockEngine(const BasisSet& _obs, const BasisSet& _dfbs) :
		obs(_obs), dfbs(_dfbs)
		{
		}

	typedef btas::RangeNd<CblasRowMajor, std::array<long, 3> > Range3d;
    	typedef btas::Tensor<double, Range3d> Tensor3d;
    	Tensor3d xyK;

    	// a DF-based builder, using coefficients of occupied MOs
    	Matrix compute_2body_fock_dfC(const Matrix& Cocc);
  	};
#endif // HAVE_DENSITY_FITTING







//----------------------------------------------------------------------------------------------------------------------
std::vector<Atom> read_geometry(const std::string& filename);
























}; /* namespace HartreeFock */


#endif /* HARTREEFOCKCLASS_H_ */

//#endif // END ADDED TO AVOID DOUBLE INCLUDE THIS HEADER FILE.
