#*********************************************************
#*                                                        *
#*      ANT.G09-2.4.1 - Makefile                          *
#*                                                        *
#**********************************************************
#*                                                        *
#*  Copyright (c) by                                      *
#*                                                        *
#*  David Jacob (1)                                       *
#*  Juan Jose Palacios (2)                                *
#*                                                        *
#* (1) Theory Department                                  *
#*     Max-Planck-Institute for Microstructure Physics    *
#*     Halle, 06120 (GERMANY)                             *
#* (2) Departamento de Fisica de la Materia Condensada    *
#*     Universidad Autonoma de Madrid                     *      
#*     28049 Madrid (SPAIN)                               *
#**********************************************************

# Makefile: A standard Makefile for libAntG
ifndef ANTGDIR
#  ANTGDIR=$(shell pwd)
  ANTGDIR=$(TOPDIR)/AntG
#  MATMULTDIR := ${CURDIR}
endif
ifndef TOPDIR
  TOPDIR=..
endif



# set complier/linker and options
#MAKE=make
FC=ifort
F90=ifort
CC=icc
AR=xiar
BIN=../bin

# MKL related libs
MKLROOT= /opt/intel/composer_xe_2015.0.090/mkl
MKLPATH= ${MKLROOT}/lib/intel64
MKLINCLUDE=${MKLROOT}/include

MODINCLUDE=$(TOPDIR)/mod

#INTELLIBS= -Wl,--start-group ${MKLPATH}/libmkl_intel_ilp64.a ${MKLPATH}/libmkl_core.a ${MKLPATH}/libmkl_intel_thread.a -Wl,--end-group
#INTELLIBS= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group
INTELDYNAMICLIBS=${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

#include make.in
# Preprocessor flags. Must be the same value as in Gaussian
#PREPROC= -DDEFMAXSHL=250000 -DDEFMAXATM=250000 -DDEFMAXTIT=100 -DDEFMAXRTE=4000
DDEFMAXSHL=250000
DDEFMAXATM=250000
DDEFMAXTIT=100
DDEFMAXRTE=4000

# Compilation flags.
#CFLAGS= -i8 -r8 -mcmodel=medium -parallel -axAVX -fPIC -shared-intel -static-libgcc -no-prec-div -fpp3 -ftz -pad -mkl -m64 -mp=nonuma -O2 -unroll # ORIGINAL ANT
#CFLAGS= -DBUILD_F90_MODULE -O3 -i8 -r8 -fpe3 -heap-arrays -mcmodel=medium -parallel -axAVX -fPIC -shared-intel -static-libgcc -no-prec-div -fpp3 -ftz -pad -mkl -m64 -mp=nonuma -unroll
#CFLAGS= -i8 -r8 -integer-size 64 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
#CFLAGS= -i8 -r8 -check all -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include # ADD -CHECK ALL FOR ERRORS CONCERNING FREE() IN DEALLOCATIONS
CFLAGS= -i8 -r8 -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include -I$(MODINCLUDE)
#CFLAGS= -i8 -r8 -I/opt/intel/composer_xe_2015.0.090/mkl/include/intel64/lp64 -I/opt/intel/composer_xe_2015.0.090/mkl/include

#CXXFLAGS+=  -I/opt/intel/composer_xe_2015.0.090/mkl/include
#FFLAGS+="-DBUILD_F90_MODULE -O3 -i8 -r8 -fpe3 -heap-arrays"
#CFLAGS =  -i8 -r8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
# Linking flags
#LFLAGS= -lpthread -lm -lc
#LDFLAGS+=   ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

#LDFLAGS+=  ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -lc -ldl -openmp -lpthread
DYNAMICLDFLAGS=${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#DYNAMICLDFLAGS=/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/libmkl_blas95_lp64.a /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/libmkl_lapack95_lp64.a -L/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
STATICLDFLAGS=${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
#STATICLDFLAGS=/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/libmkl_blas95_lp64.a /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

# use these sources:: 
#ANTSRCS= preproc.F util.f90 parameters.f90 antcommon.f90 g09Common.f90 constants.f90 numeric.f90 cluster.f90 ortho.f90 correlation.f90 OneDLead.f90 BetheLattice.f90 device.f90 ANT.f90
#ANTSRCS= echo $(ls $(ANTGDIR)/{preproc.F,util.f90,parameters.f90,antcommon.f90,g09Common.f90,constants.f90,numeric.f90 cluster.f90,ortho.f90,correlation.f90,OneDLead.f90,BetheLattice.f90,device.f90,ANT.f90})
FANTSRCS= $(ANTGDIR)/preproc.F 
F90ANTSRCS= $(ANTGDIR)/util.f90 $(ANTGDIR)/parameters.f90 $(ANTGDIR)/antcommon.f90 $(ANTGDIR)/g09Common.f90 $(ANTGDIR)/constants.f90 $(ANTGDIR)/numeric.f90 $(ANTGDIR)/cluster.f90 $(ANTGDIR)/ortho.f90 $(ANTGDIR)/correlation.f90 $(ANTGDIR)/OneDLead.f90 $(ANTGDIR)/BetheLattice.f90 $(ANTGDIR)/device.f90 $(ANTGDIR)/ANT.f90

ANTSRCS= $(FANTSRCS) $(F90ANTSRCS)

#ANTSRCS= $(ANTGDIR)/preproc.F $(ANTGDIR)/util.f90 $(ANTGDIR)/parameters.f90 $(ANTGDIR)/antcommon.f90 $(ANTGDIR)/g09Common.f90 $(ANTGDIR)/constants.f90 $(ANTGDIR)/numeric.f90 $(ANTGDIR)/cluster.f90 ortho.f90 $(ANTGDIR)/correlation.f90 $(ANTGDIR)/OneDLead.f90 $(ANTGDIR)/BetheLattice.f90 $(ANTGDIR)/device.f90 $(ANTGDIR)/ANT.f90

#SRCS= preproc.F util.f90 parameters.f90 antcommon.f90 g09Common.f90 constants.f90 lapack_blas.f numeric.f90 cluster.f90 ortho.f90 correlation.f90 OneDLead.f90 BetheLattice.f90 device.f90 ANT.f90

# use name of parent directory as library name
PARDIRNAME=$(shell dirname `pwd`)
#LIBNAME=$(shell basename $(PARDIRNAME))
LIBNAME=libAntG
EXENAME=AntG

# objects
#OBJS= $(SRCS::%.c=../obj/%.o) $(SRCS::%.f=../obj/%.o) $(SRCS::%.F=../obj/%.o) $(SRCS::%.f90=../obj/%.o) 
#ANTOBJS= $(ANTSRCS::$(ANTGDIR)/%.c=$(ANTGDIR)/%.o) $(ANTSRCS::$(ANTGDIR)/%.f=$(ANTGDIR)/%.o) $(ANTSRCS::$(ANTGDIR)/%.F=$(ANTGDIR)/%.o) $(ANTSRCS::$(ANTGDIR)/%.f90=$(ANTGDIR)/%.o)
#ANTOBJS= $(SRCS::$(ANTGDIR)/%.c=$(ANTGDIR)/%.o) $(SRCS::$(ANTGDIR)/%.f=$(ANTGDIR)/%.o) $(SRCS::$(ANTGDIR)/%.F=$(ANTGDIR)/%.o) $(SRCS::$(ANTGDIR)/%.f90=$(ANTGDIR)/%.o) 

FANTOBJS = $(subst F,o,$(FANTSRCS))
F90ANTOBJS = $(subst f90,o,$(F90ANTSRCS))
ANTOBJS = $(FANTOBJS) $(F90ANTOBJS)


# ANT.G library to be generated
ANTLIBA=$(ANTGDIR)/$(LIBNAME).a
#L502=$(BIN)/l502.exe 
ANTLIBSO=$(ANTGDIR)/$(LIBNAME).so
ANTEXE=$(ANTGDIR)/$(EXENAME).exe
ANTGF90=$(ANTGDIR)/$(EXENAME).f90



# rule for object generation
$(ANTGDIR)/%.o:: $(ANTGDIR)/%.c 
	$(CC) -c -o $@ $<
	cp $@ $(ANTGDIR)/obj
$(ANTGDIR)/%.o:: $(ANTGDIR)/%.F
	$(FC) $(CFLAGS) -I$(MKLINCLUDE) $(PREPROC) -c -o $@ $< -module $(TOPDIR)/mod
	cp $@ $(ANTGDIR)/obj
$(ANTGDIR)/%.o:: $(ANTGDIR)/%.f 
	$(FC) $(CFLAGS) -I$(MKLINCLUDE) -c -o $@ $< -module $(TOPDIR)/mod
	cp $@ $(ANTGDIR)/obj
$(ANTGDIR)/%.o:: $(ANTGDIR)/%.f90 
	$(FC) $(CFLAGS) -I$(MKLINCLUDE) -c -o $@ $< -module $(TOPDIR)/mod
	cp $@ $(ANTGDIR)/obj

# targets:: libraries and l502.exe
#target:: $(ANTLIB) $(L502)
#target:: $(ANTLIB) $(ANTGDIR)/libAntG.so $(ANTGDIR)/AntG.exe
target:: 
	echo "ANTGDIR = "$(ANTGDIR)
	echo "ANTSRCS = "$(ANTSRCS)
	echo "ANTOBJS = "$(ANTOBJS)
	make $(ANTOBJS) $(ANTLIBA) $(ANTLIBSO) $(ANTEXE)

default:: target

# rule for the library
$(ANTLIBA):: $(ANTOBJS)
	echo "ANTGDIR = "$(ANTGDIR)
	echo "ANTSRCS = "$(ANTSRCS)
	echo "ANTOBJS = "$(ANTOBJS)
	- rm  $(ANTLIBA)
	$(AR) crv $(ANTLIBA) $(ANTOBJS)
	cp $@ $(TOPDIR)/lib
#	$(AR) crv $(ANTLIB) $(ANTOBJS) $(STATICLDFLAGS)

# rule for l502.exe
#$(L502):: $(ANTOBJS)
#	$(FC) $(CFLAGS) -o $(L502) $(BIN)/ml502.o $(BIN)/l502.a $(BIN)/bdam1.o $(BIN)/caldsu.o -L$(MKLPATH) $(GAUSSLIBS) $(ANTLIB) $(INTELLIBS) $(LFLAGS) 

#libAntG.so:: $(ANTOBJS)
$(ANTLIBSO):: $(ANTOBJS)
#	$(FC) -o libAntG.so $(ANTOBJS) -L$(MKLPATH) $(INTELLIBS) $(LDFLAGS) $(ANTLIB) -shared
	$(FC) -nofor-main -o $@ $(ANTLIBA) $(DYNAMICLDFLAGS) -shared
	cp $@ $(TOPDIR)/lib

$(ANTEXE):: $(ANTOBJS)
##	$(FC) -o AntG.exe $(ANTOBJS) -L$(MKLPATH) $(INTELLIBS) $(LDFLAGS) $(ANTLIB) -shared
#	$(FC) $(CFLAGS) -c AntG.f90 -o ../obj/AntG.o
##	$(FC) -o AntG.exe AntG.o -L../lib -lAntG -L$(MKLPATH) $(INTELLIBS) $(LDFLAGS) -fPIC -openmp -lpthread
##	$(FC) -o AntG.exe AntG.o $(ANTLIB) -L$(MKLPATH) $(INTELLIBS) $(LDFLAGS) -fPIC -openmp -lpthread -static
#	$(FC) -o AntG.exe ../obj/AntG.o $(ANTLIBA) $(STATICLDFLAGS) -fPIC
	$(FC) $(CFLAGS) -c $(ANTGF90) -o $(TOPDIR)/obj/AntG.o -module $(TOPDIR)/mod
	echo "STATICLDFLAGS = "$(STATICLDFLAGS)
	$(FC) -o $@ $(TOPDIR)/obj/AntG.o $(ANTLIBA) $(STATICLDFLAGS) -fPIC
	cp $@ $(TOPDIR)/bin

#	chmod o-rwx $(L502)

# Cleaning up to compile from scratch afterwards
clean::
	@echo "Cleaning up..."
	-rm $(ANTGDIR)/*.mod
#	-rm $(TOPDIR)/*.ff90
#	-rm $(BIN)/l502.exe
	-rm $(ANTGDIR)/*.o
#	-rm $(ANTGDIR)/*.o
	-rm $(ANTGDIR)/*.a
