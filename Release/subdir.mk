################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../BowtieEntry.cpp \
../EM.cpp \
../EM_Map.cpp \
../FastaEntry.cpp \
../MarkovChain.cpp \
../Quality.cpp \
../Random_EM_Map.cpp \
../Read.cpp \
../Sift.cpp \
../main.cpp \
../test.cpp 

OBJS += \
./BowtieEntry.o \
./EM.o \
./EM_Map.o \
./FastaEntry.o \
./MarkovChain.o \
./Quality.o \
./Random_EM_Map.o \
./Read.o \
./Sift.o \
./main.o \
./test.o 

CPP_DEPS += \
./BowtieEntry.d \
./EM.d \
./EM_Map.d \
./FastaEntry.d \
./MarkovChain.d \
./Quality.d \
./Random_EM_Map.d \
./Read.d \
./Sift.d \
./main.d \
./test.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


