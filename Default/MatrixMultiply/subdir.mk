################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../MatrixMultiply/MatrixMultiplyMod.f90 


# Each subdirectory must supply rules for building sources it contributes
MatrixMultiply/%.o: ../MatrixMultiply/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) Intel(R) 64 Fortran Compiler'
	ifort -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

MatrixMultiply/MatrixMultiplyMod.o: ../MatrixMultiply/MatrixMultiplyMod.f90


