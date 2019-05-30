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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>



using namespace std;

using namespace HartreeFock;
//using HartreeFock::HartreeFockClass;

//int main(int argc, const char* argv[]) {
int main (int argc, char **argv) {

// read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
    const auto jobname = (argc > 1) ? argv[1] : "h2";
    const auto basisname = (argc > 2) ? argv[2] : "STO-3G";
    const auto fullname = (argc > 3) ? argv[3] : "h2_STO-3G_cycle0";


    int nflag = 0;
    int bflag = 0;
    int dflag = 0;
    int eflag = 0;
    char *cvalue = NULL;
    int index;
    int c;
    bool debugflag = false;
    bool erkaleflag = false;

    //char *optjobname = '\0'; // NULL DOES NOT WORK AND BREAKS STD::COUT
    //char *optbasisname = '\0'; // NULL DOES NOT WORK AND BREAKS STD::COUT
    char *optjobname = ""; // "" WORKS
    char *optbasisname = ""; // "" WORKS


    while ((c = getopt (argc, argv, "den:b:c:")) != -1)
        switch (c)
          {
          case 'd':
            dflag = 1;
            debugflag = true;
            break;
          case 'e':
            eflag = 1;
            erkaleflag = true;
            break;
          case 'n':
            nflag = 1;
            optjobname = optarg;
            break;
          case 'b':
            bflag = 1;
            optbasisname = optarg;
            break;
          case 'c':
            cvalue = optarg;
            break;
          case '?':
            if (optopt == 'c')
              fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
              fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
              fprintf (stderr,
                       "Unknown option character `\\x%x'.\n",
                       optopt);
            return 1;
          default:
            abort ();
          }


    printf ("nflag = %d, bflag = %d, dflag = %d, eflag = %d, cvalue = %s\n",
              nflag, bflag, dflag,eflag, cvalue);

    for (index = optind; index < argc; index++)
      printf ("Non-option argument %s\n", argv[index]);


    cout <<"jobname: "<< jobname <<"\n";
    cout <<"basisname: "<< basisname <<"\n";
    cout <<"optjobname: "<< optjobname <<"\n";
    cout <<"optbasisname: "<< optbasisname <<"\n";

    cout <<"fullname: "<< fullname <<"\n";


	HartreeFockClass hf;
	//hf.run("h2o.xyz","aug-cc-pVDZ");

	//hf.runant("h2","STO-3G");
	hf.runant(jobname,basisname,fullname,debugflag,erkaleflag);

	cout <<"Exited hf.runant(jobname,basisname,fullname,debugflag,erkaleflag) Successfully\n";
	cout <<"Exited hf.runant("<<jobname<<","<<basisname<<","<<fullname<<","<<(debugflag ? "true" : "false")<<","<<(erkaleflag ? "true" : "false")<<") Successfully\n";

		size_t N = 3;
		vector<double> A(N * N), B(N * N), C(N * N);

//		HartreeFock::HartreeFockClass hf;

	if(debugflag == true){
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


		//hf.mat_mult(A, sz, sz, B, sz, sz, C, sz, sz);


		//print A, B, C on the standard output
		hf.print_matrix(A, sz, sz);
		hf.print_matrix(B, sz, sz);
		hf.print_matrix(C, sz, sz);
	}


	//return 0;


	return 0;
}
