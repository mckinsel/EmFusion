################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../BowtieEntry.cpp \
../EM.cpp \
../EM_Map.cpp \
../Quality.cpp \
../Read.cpp 

OBJS += \
./BowtieEntry.o \
./EM.o \
./EM_Map.o \
./Quality.o \
./Read.o 

CPP_DEPS += \
./BowtieEntry.d \
./EM.d \
./EM_Map.d \
./Quality.d \
./Read.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


