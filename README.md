NEGFY is a Nonequilibrium Green's Function HF/DFT Quantum Electronics code based in the NEGF formalism in its matrix form on a basis of Gaussian Orbitals.

NEGFY merges three Quantum Physics codes, exploiting their particular capabilities:

1. ANT.G (https://www.simuneatomistics.com/ant/): developed by JJ Palacios and coworkers, computes the electronic structure with open boundary conditions in equilibrium and non-equilibrium situations (i.e., zero and finite bias voltage). It follows an onion shell structure, modelling the far electrodes with a parametrized Bethe lattice while computing at  the density functional theory level, as implemented in GAUSSIAN03/09, the central scattering region. Used through NEGFY, ANT is capable to work without calling GAUSSIAN.

2. Libint (https://github.com/evaleev/libint): Evgeny Valeev's library for the evaluation of molecular integrals of many-body operators over Gaussian functions, provides the matrix form of al One-Body and Two-Body operators at the Hartree-Fock level.

3. ERKALE -- HF/DFT from Helgaker (https://github.com/susilehtola/erkale): Susi Lehtola's DFT code provides the matrix form of al One-Body and Two-Body operators at the DFT level, by interfacing both Libint and the DFT Functional library LibXC.

· e-mail - carlos.salgado@uam.es

Prerequisites:

· Libint (https://github.com/evaleev/libint).

· LibXC (https://tddft.org/programs/libxc): a library of exchange-correlation functionals for density-functional theory. The aim is to provide a portable, well tested and reliable set of exchange and correlation functionals that can be used by all the ETSF codes and also other codes. In Libxc you can find different types of functionals: LDA, GGA, hybrids, and mGGA. These functionals depend on local information, in the sense that the value of the potential at a given point depends only on the values of the density – and the gradient of the density and the kinetic energy density, for the GGA and mGGA cases – at a given point. It can calculate the functional itself and its derivative; for most functionals, higher-order derivatives are available. Libxc is written in C and has Fortran and Python bindings.

· Eigen (http://eigen.tuxfamily.org): a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.

· Armadillo (http://arma.sourceforge.net): a high quality linear algebra library (matrix maths) for the C++ language, aiming towards a good balance between speed and ease of use.

· A C++ and Fortran compiler, preferably the Intel Compiler. GCC compiler is also tested.

· CMake v3.4.3 or higher (https://cmake.org/): CMake is an open-source, cross-platform family of tools designed to build, test and package software.

· HDF5 (https://www.hdfgroup.org/solutions/hdf5/): High-performance data management and storage suite.

Under construction...
