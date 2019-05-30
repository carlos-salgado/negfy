# Makefile: A standard Makefile for HartreeFockClass.cpp and others.
ifndef TOPDIR
  TOPDIR=..
endif
ifndef HFDIR
  HFDIR=$(TOPDIR)/HartreeFock
endif

LIBINTROOT=/home/carlos/Python-Projects/DFT-CODES-PYTHON/libint-2.2-spher-shared-intel-derivs
LIBINTBUILDDIR=$(LIBINTROOT)/libint-build
LIBINTMASTERDIR=$(LIBINTROOT)/libint-master
-include $(LIBINTBUILDDIR)/src/bin/MakeVars
-include $(LIBINTBUILDDIR)/src/lib/libint/MakeVars.features
CPP = icpc
CXX = icpc
AR = ar
CPPFLAGS += -I$(LIBINTBUILDDIR)/include -I$(LIBINTBUILDDIR)/include/libint2 -I$(LIBINTBUILDDIR)/src/lib/libint -DSRCDATADIR=\"$(LIBINTMASTERDIR)/lib/basis\" 
LDFLAGS  =
FORTLIB  = -L$(HFDIR)/../lib -lMatrixMultiplyMod
FORTSTATICLIB  = $(HFDIR)/../lib/libMatrixMultiplyMod.a


COMPILER_LIB = $(LIBINTBUILDDIR)/src/bin/libint/libINT.a
#COMPUTE_LIB = -lint2
#COMPUTE_LIB = '-L$(LIBINTBUILDDIR)/lib -lint2'
COMPUTE_LIB = $(LIBINTBUILDDIR)/lib/libint2.so

LINKTOCOMPUTE_LIB = '-L$(LIBINTBUILDDIR)/lib/libint2.so -lint2'
vpath %.a $(LIBINTBUILDDIR)/lib:$(LIBINTBUILDDIR)/lib/.libs

OBJSUF = o
DEPSUF = d
CXXDEPENDSUF = none
CXXDEPENDFLAGS = -M

default:: all
#default:: target

#all:: HartreeFockClass.o libHartreeFock.a libHartreeFock.so HartreeFock.exe $(TEST2) 
#all:: $(TEST2) libHartreeFock.a libHartreeFock.so HartreeFock.exe
#all::
#	cd $(HFDIR)
#	echo "PWD ="$(shell pwd)
#	echo "HFDIR ="$(HFDIR)
#	echo "make target inside HartreeFock"
#	make target
#	cd $(HFDIR)/..

#target::
#	cd $(HFDIR)
#	echo "PWD ="$(shell pwd)
#	echo "HFDIR ="$(HFDIR)
#	echo "make target inside HartreeFock"
#	make things
#	cd $(HFDIR)/..

TEST2 = $(HFDIR)/HartreeFockClass
CXXTEST2SRC = $(TEST2).cpp
CXXTEST2OBJ = $(CXXTEST2SRC:%.cpp=%.$(OBJSUF))
CXXTEST2DEP = $(CXXTEST2SRC:%.cpp=%.$(DEPSUF))

HFEXE = $(HFDIR)/HartreeFock.exe
CXXHFEXESRC = $(HFDIR)/main.cpp
#CXXHFEXEOBJ = $(CXXTEST2SRC:%.cpp=%.$(OBJSUF)) main.o
CXXHFEXEOBJ = $(CXXHFEXESRC:%.cpp=%.$(OBJSUF))
#CXXHFEXEDEP = $(CXXTEST2SRC:%.cpp=%.$(DEPSUF))
CXXHFEXEDEP = 

HFLIBA = $(HFDIR)/libHartreeFock.a
#CXXTEST2SRC = $(TEST2).cpp
#CXXTEST2OBJ = $(CXXTEST2SRC:%.cpp=%.$(OBJSUF))
#CXXTEST2DEP = $(CXXTEST2SRC:%.cpp=%.$(DEPSUF))

HFLIBSO = $(HFDIR)/libHartreeFock.so
#CXXTEST2SRC = $(TEST2).cpp
#CXXTEST2OBJ = $(CXXTEST2SRC:%.cpp=%.$(OBJSUF))
#CXXTEST2DEP = $(CXXTEST2SRC:%.cpp=%.$(DEPSUF))

#target:: $(TEST2) $(HFDIR)/libHartreeFock.a $(HFDIR)/libHartreeFock.so $(HFDIR)/HartreeFock.exe
#target:: $(TEST2) $(HFDIR)/libHartreeFock.a $(HFDIR)/libHartreeFock.so $(HFEXE)
#target:: $(HFDIR)/libHartreeFock.a $(HFDIR)/libHartreeFock.so $(HFEXE)
target:: 
	cd $(HFDIR)
	make LIBS
	cd $(TOPDIR)

LIBS:: $(HFLIBA) $(HFLIBSO) $(HFEXE)

#things:: $(TEST2) libHartreeFock.a libHartreeFock.so HartreeFock.exe


$(TEST2):: $(CXXTEST2OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS) $(FORTLIB) -lpthread
	cp $@ $(TOPDIR)/bin

#$(HFLIBA):: $(CXXTEST2OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
$(HFLIBA):: $(CXXTEST2OBJ)
	$(AR) rc $@ $^ 
	cp $@ $(TOPDIR)/lib

$(HFLIBSO):: $(CXXTEST2OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD)  -o $@ $(LDFLAGS) $^ $(SYSLIBS) $(FORTLIB) -lpthread -shared
	cp $@ $(TOPDIR)/lib

$(HFEXE):: $(CXXHFEXEOBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS) $(FORTLIB) $(HFLIBA) -lpthread
	cp $@ $(TOPDIR)/bin


#$(HFDIR)/libHartreeFock.a::
##	$(AR) rc $@ $(LDFLAGS) $^ $(FORTSTATICLIB)
#	$(LD)  -o $@ $(LDFLAGS) $^ $(SYSLIBS) $(FORTSTATICLIB) -static-intel -nostartfiles

#$(HFDIR)/libHartreeFock.so::
#	$(LD)  -o $@ $(LDFLAGS) $^ $(SYSLIBS) $(FORTLIB) -lpthread -shared

#$(HFDIR)/HartreeFock.exe::
#	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS) $(FORTLIB) -lpthread -fPIC


## Source files for timer and tester are to be compiled using CXXGEN
#$(TEST2):: CXX=$(CXXGEN)
#$(TEST2):: CXXFLAGS=$(CXXGENFLAGS)
#$(TEST2):: LD=$(CXXGEN)

#$(HFDIR)/libHartreeFock.a $(HFDIR)/libHartreeFock.so $(HFDIR)/HartreeFock.exe:: CXX=$(CXXGEN)
#$(HFDIR)/libHartreeFock.a $(HFDIR)/libHartreeFock.so $(HFDIR)/HartreeFock.exe:: CXXFLAGS=$(CXXGENFLAGS)
#$(HFDIR)/libHartreeFock.a $(HFDIR)/libHartreeFock.so $(HFDIR)/HartreeFock.exe:: LD=$(CXXGEN)

# Source files for timer and tester are to be compiled using CXXGEN
$(TEST2) $(HFLIBA) $(HFLIBSO) $(HFEXE):: CXX=$(CXXGEN)
$(TEST2) $(HFLIBA) $(HFLIBSO) $(HFEXE):: CXXFLAGS=$(CXXGENFLAGS)
$(TEST2) $(HFLIBA) $(HFLIBSO) $(HFEXE):: LD=$(CXXGEN)


clean::
	-rm -rf $(TEST1) $(TEST2) *.o *.d *.a *.so *.exe

distclean:: realclean
	-rm -rf $(LIBINTBUILDDIR)/include/libint2/boost

realclean:: clean

targetclean:: clean

$(LIBINTBUILDDIR)/include/libint2/boost/preprocessor.hpp: $(LIBINTBUILDDIR)/external/boost.tar.gz
	gunzip -c $(LIBINTMASTERDIR)/external/boost.tar.gz | tar -xf - -C $(LIBINTBUILDDIR)/include/libint2

depend:: $(CXXTEST2DEP) $(CXXHFEXEDEP)

ifneq ($(CXXDEPENDSUF),none)
%.d:: %.cpp $(LIBINTBUILDDIR)/include/libint2/boost/preprocessor.hpp
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< > /dev/null
	sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' < $(*F).$(CXXDEPENDSUF) > $(@F)
	/bin/rm -f $(*F).$(CXXDEPENDSUF)
else
%.d:: %.cpp $(LIBINTBUILDDIR)/include/libint2/boost/preprocessor.hpp
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< | sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' > $(@F)
	cp $@ $(TOPDIR)/obj
endif

#-include $(CXXTEST1DEP)
-include $(CXXTEST2DEP) $(CXXHFEXEDEP)

#$(HFDIR)/HartreeFockClass.o::
##	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< | sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' > $(@F)
#	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< | find $(HFDIR) -type f -exec sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' < $(*F).$(CXXDEPENDSUF) > $(@F)



