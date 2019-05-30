//============================================================================
// Name        : 20170321-hartree-fock-antg-cpp-fortran.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++
//============================================================================

#include "HartreeFock/HartreeFockClass.h"
#include <iostream>
#include <vector>
#include <random>
#include <ctime>


using namespace std;

using namespace HartreeFock;

int main(void) {
	cout <<"hello world";

	size_t N = 3;
	vector<double> A(N * N), B(N * N), C(N * N);

	HartreeFock::HartreeFockClass rnd;
	rnd.run("h2o.xyz","aug-cc-pVDZ");
//        const auto filename = "h2o.xyz";
//        const auto basisname = "aug-cc-pVDZ";


	//Fill A and B with random numbers in the range [0,1]:
	rnd.rnd_fill(A, 0.0, 1.0);
	rnd.rnd_fill(B, 0.0, 1.0);

	//Multiply matrices A and B, save the result in C
	int sz = N;
	//rnd.matrix_multiply(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
	//rnd.mat_mult(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
	//rnd.mat_mult(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
	//rnd.mat_mult(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
	//rnd.mat_mult(&A[0], &sz, &sz, &B[0], &sz, &sz, &C[0], &sz, &sz);
	rnd.mat_mult(A, sz, sz, B, sz, sz, C, sz, sz);


	//print A, B, C on the standard output
	rnd.print_matrix(A, sz, sz);
	rnd.print_matrix(B, sz, sz);
	rnd.print_matrix(C, sz, sz);

	return 0;
}
