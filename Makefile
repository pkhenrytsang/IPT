# C++ compiler
CXX = icpc
#CXX = g++-8
CC = gcc
MPICXX = mpiicpc

prefix = ~/bin

# CUDA compiler

# CUDA Libraries
CUDAPATH = /usr/local/cuda/
#nvcc

gpuprog = IPT-gpu
cpuprog = IPT-cpu
mcprog = IPT-mc
main = IPT
mainmc = IPTMC


# source path
SP = src
#object files path
OP = obj
# executable path
RP = bin

#mkdir
MKDIR_P = mkdir -vp


# --- Flags (GNU) --- #


#Flags for GNU C++
#CXXFLAGS = -fopenmp -march=native -O2 $(INC) -std=c++14

#Flags for CUDA (GNU)
#NVCCFLAGS = -arch=sm_75 -ccbin $(CXX) -Xcompiler -fopenmp -O2 -std=c++14


# --- Flags (Intel) --- #

#Flags for Intel C++
CXXFLAGS = -Wall -qopenmp -xHost -O2 $(INC) -std=c++14

MPICXXFLAGS = -Wall -xHost -O2 -std=c++14

#Flags for CUDA (Intel)
NVCCFLAGS = -arch=sm_75 -ccbin $(CXX) -Xcompiler -qopenmp -O2 -std=c++14

# --- Libs --- #

#LIBS for Intel C++
LIBS = -mkl -lgsl -lgslcblas

#LIBS for GNU C++
#LIBS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lgsl -lgslcblas

#CUDA-related
CUDALIBS =  -L$(CUDAPATH)lib64 -lcuda -lcudart
INC = -I$(CUDAPATH)include

all : directories gpuprogram cpuprogram mcprogram utilities

gpu: directories gpuprogram

cpu: directories cpuprogram

mc: directories mcprogram

utilities : directories analytic_continuation

gpuprogram : $(OP)/$(main).o $(OP)/SIAM.o $(OP)/Grid.o $(OP)/Params.o $(OP)/routines.o $(OP)/dinterpl.o $(OP)/SIAM_GPU.o $(OP)/tail.o
	$(CXX) $(CXXFLAGS) -o $(RP)/$(gpuprog) $(OP)/$(main).o $(OP)/SIAM.o $(OP)/Grid.o $(OP)/Params.o $(OP)/routines.o $(OP)/dinterpl.o $(OP)/SIAM_GPU.o $(OP)/tail.o $(LIBS) $(CUDALIBS)
	
cpuprogram : $(OP)/$(main).o $(OP)/SIAM.o $(OP)/Grid.o $(OP)/Params.o $(OP)/routines.o $(OP)/dinterpl.o $(OP)/SIAM_CPU.o $(OP)/tail.o
	$(CXX) $(CXXFLAGS) -o $(RP)/$(cpuprog) $(OP)/$(main).o $(OP)/SIAM.o $(OP)/Grid.o $(OP)/Params.o $(OP)/routines.o $(OP)/dinterpl.o $(OP)/SIAM_CPU.o $(OP)/tail.o $(LIBS)
	
mcprogram : $(OP)/$(mainmc).o $(OP)/SIAMMC.o $(OP)/GridMC.o $(OP)/Params.o $(OP)/routines.o $(OP)/dinterpl.o $(OP)/tail.o
	$(MPICXX) $(MPICXXFLAGS) -o $(RP)/$(mcprog) $(OP)/$(mainmc).o $(OP)/SIAMMC.o $(OP)/GridMC.o $(OP)/Params.o $(OP)/routines.o $(OP)/dinterpl.o $(OP)/tail.o $(LIBS)

analytic_continuation : $(OP)/acond.o $(OP)/routines.o
	$(CXX) $(CXXFLAGS) -o $(RP)/acond $(OP)/acond.o $(OP)/routines.o $(LIBS)

directories : $(OP) $(RP)

#Create directories
$(OP) :
	$(MKDIR_P) $(OP)

$(RP) :
	$(MKDIR_P) $(RP)
	
# main program
$(OP)/$(main).o : $(SP)/$(main).cpp $(SP)/SIAM.h $(SP)/Grid.h
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/$(main).cpp
	
# main program
$(OP)/$(mainmc).o : $(SP)/$(mainmc).cpp $(SP)/SIAMMC.h $(SP)/GridMC.h
	$(MPICXX) $(MPICXXFLAGS) -c -o $@ $(SP)/$(mainmc).cpp

# SIAM
$(OP)/SIAM.o : $(SP)/SIAM.cpp $(SP)/SIAM.h $(SP)/Grid.h $(SP)/routines.h
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/SIAM.cpp
	
# SIAM
$(OP)/SIAMMC.o : $(SP)/SIAMMC.cpp $(SP)/SIAMMC.h $(SP)/GridMC.h $(SP)/routines.h
	$(MPICXX) $(MPICXXFLAGS) -c -o $@ $(SP)/SIAMMC.cpp
	
# SIAM
$(OP)/SIAM_CPU.o : $(SP)/SIAM.cpu.cpp $(SP)/SIAM.h $(SP)/Grid.h $(SP)/routines.h
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/SIAM.cpu.cpp

# cuSIAM (GPU)
$(OP)/SIAM_GPU.o : $(SP)/SIAM.cu $(SP)/SIAM.h
	nvcc $(NVCCFLAGS) -c -o $@ $(SP)/SIAM.cu

# Grid
$(OP)/Grid.o : $(SP)/Grid.cpp $(SP)/Grid.h
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/Grid.cpp
	
# Grid
$(OP)/GridMC.o : $(SP)/GridMC.cpp $(SP)/GridMC.h
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/GridMC.cpp

# Input class used for reading files with parameters
$(OP)/Params.o : $(SP)/Params.cpp $(SP)/Params.h
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/Params.cpp

# Interpolation
$(OP)/dinterpl.o : $(SP)/dinterpl.cpp $(SP)/dinterpl.h
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/dinterpl.cpp

# contains some constants and useful numerical routines
$(OP)/routines.o : $(SP)/routines.cpp $(SP)/routines.h 
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/routines.cpp

# Analytic continuation to Matsubara axis
$(OP)/acond.o : $(SP)/acond.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/acond.cpp


$(OP)/tail.o : $(SP)/tail.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $(SP)/tail.cpp
	
# clean all object and exec files
clean :
	rm -vf $(RP)/* $(OP)/*.o
	rm -vr $(RP) $(OP)

install :
	cp -v $(RP)/* $(prefix)/.
