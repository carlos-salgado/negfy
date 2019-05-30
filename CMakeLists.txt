#cmake_minimum_required(VERSION 2.8.9)
cmake_minimum_required (VERSION 2.6)
project(directory_test)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix "." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the NEGFY directory."
      FORCE)

# ADDED TO LINK WITH "ar crv" instead of "ar cq", WHICH FAILS TO LINK EXECUTABLES.
# SEEMS TO BE NOT NECESSARY.
#SET(CMAKE_C_ARCHIVE_CREATE "echo \"Top directory\"; <CMAKE_AR> crv <TARGET> <LINK_FLAGS> <OBJECTS>")
#SET(CMAKE_CXX_ARCHIVE_CREATE "echo \"Top directory\"; <CMAKE_AR> crv <TARGET> <LINK_FLAGS> <OBJECTS>")
#SET(CMAKE_Fortran_ARCHIVE_CREATE "echo \"Top directory\"; <CMAKE_AR> crv <TARGET> <LINK_FLAGS> <OBJECTS>")

set(MODINCLUDE "${default_prefix}/mod")
#set(MODINCLUDE "${CMAKE_INSTALL_PREFIX}/mod")
set(INCLUDE "${default_prefix}/include")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/usr/include -I/usr/local/include -I/usr/include/c++/5 -I${MKLROOT}/include -I${INCLUDE}")

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-i8 -r8 -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include -O3 -gen-interfaces -module ${MODINCLUDE} ")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-i8 -r8 -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include -O0 -g -gen-interfaces -module ${MODINCLUDE}")
#  set (CMAKE_Fortran_FLAGS_RELEASE " -i8 -r8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -module= ${MODINCLUDE} -O3 -gen-interfaces")
#  set (CMAKE_Fortran_FLAGS_DEBUG   " -i8 -r8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -module= ${MODINCLUDE} -O0 -g -gen-interfaces")

elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32 -module= ${MODINCLUDE}")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32 -module= ${MODINCLUDE}")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -module= ${MODINCLUDE}")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -module= ${MODINCLUDE}")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

set(MKLROOT $ENV{MKLROOT})






set(AR "ar") # COMMENTED ON 2018-04-18 TO AVOID LINKING PROBLEMS WHEN LINKING negfy.f90
#set(ARFLAGS "r")
set(ARFLAGS " --target=elf64-x86-64")
# /usr/bin/ar: objetivos admitidos: elf64-x86-64 elf32-i386 elf32-iamcu elf32-x86-64 a.out-i386-linux pei-i386 pei-x86-64 elf64-l1om elf64-k1om elf64-little elf64-big elf32-little elf32-big pe-x86-64 pe-bigobj-x86-64 pe-i386 plugin srec symbolsrec verilog tekhex binary ihex
set(LD "/usr/bin/ld -m elf_x86_64") # ORIGINAL CODE COMMENTED ON 2018-04-18
#set(LD "/usr/bin/ld -m elf64-x86-64")







set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lcilkrts -lifcore -lifport -limf -lintlc -lirng -lsvml -lm -ldl")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_cdft_core.a ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lpthread -lm -ldl")

set (CMAKE_Fortran_EXE_LINKER_FLAGS "${CMAKE_Fortran_EXE_LINKER_FLAGS} ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lcilkrts -lifcore -lifport -limf -lintlc -lirng -lsvml -lm -ldl")

#set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl")



set(LIBINTROOT "/home/carlos/Python-Projects/DFT-CODES-PYTHON/libint-2.2-spher-shared-intel-derivs")
set(LIBINTBUILDDIR "${LIBINTROOT}/libint-build")
set(LIBINTMASTERDIR "${LIBINTROOT}/libint-master")
set(HFDIR "${CMAKE_CURRENT_SOURCE_DIR}/HartreeFock")
set(SRCDATADIR "${LIBINTMASTERDIR}/lib/basis")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11 -I/usr/local/include/eigen3 -O3  -DHAVE_CONFIG_H  -I${HFDIR} -I${LIBINTMASTERDIR}/src/bin -I${LIBINTMASTERDIR}/include -I${LIBINTBUILDDIR}/include -I${LIBINTBUILDDIR}/src/bin -fPIC -std=c++11 -I/usr/local/include/eigen3 -I${LIBINTBUILDDIR}/include -I${LIBINTBUILDDIR}/include/libint2 -I${LIBINTBUILDDIR}/src/lib/libint -DSRCDATADIR=\\\"/home/carlos/Python-Projects/DFT-CODES-PYTHON/libint-2.2-spher-shared-intel-derivs/libint-master/tests/hartree-fock/../../lib/basis\\\" ")







# Add the armadillo include directory
INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})
LINK_DIRECTORIES(${ARMADILLO_LIBRARY_DIRS})
# Add the GSL include directory
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
LINK_DIRECTORIES(${GSL_LIBRARY_DIRS})
# Add libxc include directory
INCLUDE_DIRECTORIES(${LIBXC_INCLUDE_DIRS})
LINK_DIRECTORIES(${LIBXC_LIBRARY_DIRS})
# Add libint include directory
INCLUDE_DIRECTORIES(${LIBINT_INCLUDE_DIRS})
LINK_DIRECTORIES(${LIBINT_LIBRARY_DIRS})
# Add HDF5 include directory
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
LINK_DIRECTORIES(${HDF5_LIBRARY_DIRS})




set(INCLUDEDIRS "${LIBINTBUILDDIR}/src/bin" "${LIBINTBUILDDIR}/include" "${LIBINTMASTERDIR}/include" "${LIBINTMASTERDIR}/src/bin" "/usr/local/include/eigen3" "/usr/include/eigen3" "/home/carlos/Python-Projects/DFT-CODES-PYTHON/libint-2.2-spher-shared-intel-derivs/libint-build/include" "/home/carlos/Python-Projects/DFT-CODES-PYTHON/libint-2.2-spher-shared-intel-derivs/libint-build/include/libint2" "/home/carlos/Python-Projects/DFT-CODES-PYTHON/libint-2.2-spher-shared-intel-derivs/libint-build/src/lib/libint")
#Bring the headers, such as Student.h into the project
#include_directories(include)
include_directories(${INCLUDEDIRS})

 
#However, the file(GLOB...) allows for wildcard additions:
#file(GLOB NEGFYSOURCES "MatrixMultiply/*.c" "MatrixMultiply/*.f" "MatrixMultiply/*.f90")



#get_filename_component (default_prefix ".")
# EXECUTABLES
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
# STATIC LIBRARIES
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
# SHARED LIBRARIES
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")




# Find LAPACK
find_package(LAPACK REQUIRED)
message("LAPACK headers are in ${LAPACK_INCLUDE_DIRS}")
message("LAPACK library is ${LAPACK_LIBRARIES}")
message("")

# Find Armadillo
find_package(Armadillo REQUIRED)
message("ARMADILLO headers are in ${ARMADILLO_INCLUDE_DIRS}")
message("ARMADILLO library is ${ARMADILLO_LIBRARIES}")
message("")

# Find GSL
find_package(GSL REQUIRED)
message("GSL headers are in ${GSL_INCLUDE_DIRS}")
message("GSL library is ${GSL_LIBRARIES}")
message("")

# Find libxc
find_package(libxc REQUIRED)
message("libxc headers are in ${LIBXC_INCLUDE_DIRS}")
message("libxc linker flags are ${LIBXC_LIBRARIES}")
message("")

## Find libint
#find_package(libint REQUIRED)
#message("libint headers are in ${LIBINT_INCLUDE_DIRS}")
#message("libint linker flags are ${LIBINT_LIBRARIES}")
#message("")

# Find HDF5
find_package(HDF5 REQUIRED)
message("HDF5 headers are in ${HDF5_INCLUDE_DIRS}")
message("HDF5 linker flags are ${HDF5_LIBRARIES}")
message("")





add_subdirectory(MatrixMultiply)
add_subdirectory(Ant)
#add_library(negfy ${other_srcs})
add_library(negfy "negfy.f90")
target_link_libraries(negfy MatrixMultiplyMod Ant AntSt)

set(EXECUTABLENAME "negfy.exe")
add_executable(${EXECUTABLENAME} "negfy.f90")
target_link_libraries("negfy.exe" MatrixMultiplyMod Ant AntSt)


# Link libraries
if(BUILD_SHARED_LIBS)
else()
 if(UNIX AND NOT APPLE) 
  target_link_libraries(${EXECUTABLENAME} -lrt)
 endif()
 target_link_libraries(${EXECUTABLENAME} ${ARMADILLO_LIBRARIES})
 target_link_libraries(${EXECUTABLENAME} ${GSL_LIBRARIES})
 target_link_libraries(${EXECUTABLENAME} ${LAPACK_LIBRARIES})
 target_link_libraries(${EXECUTABLENAME} ${LIBXC_LIBRARIES})
 target_link_libraries(${EXECUTABLENAME} ${LIBINT_LIBRARIES})
 target_link_libraries(${EXECUTABLENAME} ${HDF5_LIBRARIES})
endif()





# ADD libint2.so and libAnt.so TO THE CMAKE_EXE_LINKER_FLAGS.
set(NEGFYLIBDIR "${CMAKE_BINARY_DIR}/lib")
set(LIBINTROOT "/home/carlos/Python-Projects/DFT-CODES-PYTHON/libint-2.2-spher-shared-intel-derivs")
set(LIBINTBUILDDIR "${LIBINTROOT}/libint-build")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lcilkrts -lifcore -lifport -limf -lintlc -lirng -lsvml -lm -ldl ${LIBINTBUILDDIR}/src/bin/libint/libINT.a ${LIBINTBUILDDIR}/lib/libint2.so ${NEGFYLIBDIR}/libAnt.so -lm -lmpfr -lgmpxx -lgmp -lpthread")

add_subdirectory(HartreeFock)
set(HFEXECUTABLENAME "negfyHF.exe")
add_executable(${HFEXECUTABLENAME} "main.cpp")
target_link_libraries("negfyHF.exe" Ant HartreeFock) # MatrixMultiplyMod)





# Link libraries
#if(BUILD_SHARED_LIBS)
#else()
# if(UNIX AND NOT APPLE) 
  target_link_libraries(${HFEXECUTABLENAME} -lrt)
# endif()
 target_link_libraries(${HFEXECUTABLENAME} ${ARMADILLO_LIBRARIES})
 target_link_libraries(${HFEXECUTABLENAME} ${GSL_LIBRARIES})
 target_link_libraries(${HFEXECUTABLENAME} ${LAPACK_LIBRARIES})
 target_link_libraries(${HFEXECUTABLENAME} ${LIBXC_LIBRARIES})
 target_link_libraries(${HFEXECUTABLENAME} ${LIBINT_LIBRARIES})
 target_link_libraries(${HFEXECUTABLENAME} ${HDF5_LIBRARIES})
#endif()



#Generate the static library from the sources
#add_library(negfy STATIC ${MATMULTSOURCES})
#Generate the static library from the sources
#add_library(negfySh SHARED ${MATMULTSOURCES})
#Generate the executable from the sources
#add_executable ("negfy.exe" ${MATMULTSOURCES})
 
#Set the location for library installation -- i.e., /usr/lib in this case
# not really necessary in this example. Use "sudo make install" to apply
install(TARGETS negfy DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
