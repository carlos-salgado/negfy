//============================================================================
// Name        : 20170321-hartree-fock-antg-cpp-fortran.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++
//============================================================================

#include "HartreeFockClass.h"
#include <iostream>
#include <vector>
#include <random>
#include <ctime>


using namespace std;

using namespace HartreeFock;
//using HartreeFock::HartreeFockClass;

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
