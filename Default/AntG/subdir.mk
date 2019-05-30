################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../AntG/ANT.f90 \
../AntG/AntG.f90 \
../AntG/AntG2.f90 \
../AntG/BetheLattice.f90 \
../AntG/BetheLattice_no1DBL.f90 \
../AntG/OneDLead.f90 \
../AntG/antcommon.f90 \
../AntG/cluster.f90 \
../AntG/cluster_no1DBL.f90 \
../AntG/constants.f90 \
../AntG/correlation.f90 \
../AntG/device.f90 \
../AntG/fortranc.f90 \
../AntG/g09Common.f90 \
../AntG/numeric.f90 \
../AntG/numeric_no1DBL.f90 \
../AntG/ortho.f90 \
../AntG/parameters.f90 \
../AntG/util.f90 

C_SRCS += \
../AntG/system.c 

C_DEPS += \
./AntG/system.d 


# Each subdirectory must supply rules for building sources it contributes
AntG/%.o: ../AntG/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) Intel(R) 64 Fortran Compiler'
	ifort -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

AntG/ANT.o: ../AntG/ANT.f90 mod/antcommon.o mod/constants.o mod/device.o mod/g09common.o AntG/parameters.o mod/preproc.o mod/util.o

AntG/AntG.o: ../AntG/AntG.f90 mod/antmod.o

AntG/AntG2.o: ../AntG/AntG2.f90

AntG/BetheLattice.o: ../AntG/BetheLattice.f90 mod/antcommon.o mod/cluster.o mod/cluster.o mod/constants.o mod/g09common.o mod/numeric.o mod/numeric.o AntG/parameters.o

AntG/BetheLattice_no1DBL.o: ../AntG/BetheLattice_no1DBL.f90 mod/antcommon.o mod/cluster.o mod/cluster.o mod/constants.o mod/g09common.o mod/numeric.o mod/numeric.o AntG/parameters.o

AntG/OneDLead.o: ../AntG/OneDLead.f90 mod/constants.o mod/g09common.o mod/numeric.o mod/numeric.o AntG/parameters.o

AntG/antcommon.o: ../AntG/antcommon.f90

AntG/cluster.o: ../AntG/cluster.f90 mod/antcommon.o mod/g09common.o AntG/parameters.o mod/preproc.o

AntG/cluster_no1DBL.o: ../AntG/cluster_no1DBL.f90 mod/antcommon.o mod/g09common.o AntG/parameters.o mod/preproc.o

AntG/constants.o: ../AntG/constants.f90

AntG/correlation.o: ../AntG/correlation.f90 mod/constants.o mod/numeric.o mod/numeric.o AntG/parameters.o mod/util.o

AntG/device.o: ../AntG/device.f90 mod/bethelattice.o mod/bethelattice.o mod/onedlead.o mod/antcommon.o mod/cluster.o mod/cluster.o mod/constants.o mod/correlation.o mod/g09common.o mod/numeric.o mod/numeric.o mod/orthogonalization.o AntG/parameters.o mod/preproc.o mod/util.o

AntG/fortranc.o: ../AntG/fortranc.f90

AntG/g09Common.o: ../AntG/g09Common.f90 AntG/parameters.o mod/preproc.o

AntG/numeric.o: ../AntG/numeric.f90 mod/antcommon.o mod/constants.o

AntG/numeric_no1DBL.o: ../AntG/numeric_no1DBL.f90 mod/antcommon.o mod/constants.o

AntG/ortho.o: ../AntG/ortho.f90 mod/constants.o mod/numeric.o mod/numeric.o mod/util.o

AntG/parameters.o: ../AntG/parameters.f90 mod/preproc.o mod/util.o

AntG/%.o: ../AntG/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

AntG/util.o: ../AntG/util.f90


