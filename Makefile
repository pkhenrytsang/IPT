# C++ compiler
Cpp = icpc 

main = cudaIPT

#mpiCC
# source path
SP = src
#object files path
OP = obj
# executable path
RP = bin

# ----------- GPU+OPENMP -----------------------------------# 

mpiCC = icpc 
FLAGS =   -qopenmp -D_OMP -xHost -O2 $(INC) 
CUDAFLAGS = -arch=sm_75 -ccbin icpc -Xcompiler -qopenmp -O2 -std=c++14

# -----------  HYBRID ------------------------------------#

#mpiCC = /opt/openmpi-1.6.3/bin/mpicxx
#FLAGS =  -D_MPI -D_OMP -openmp -static-intel #-fast

#---------------------------------------------------------#

LIBS = -mkl -lgsl -lgslcblas # use this if needed

CUDALIBS =  -L/usr/local/cuda/lib64 -lcuda -lcudart

INC = -I/usr/local/cuda/include

all : $(OP)/$(main).o $(OP)/SIAM.o $(OP)/Grid.o $(OP)/log.o $(OP)/Params.o $(OP)/routines.o $(OP)/dinterpl.o $(OP)/SIAM_GPU.o
	$(mpiCC) $(FLAGS) -o $(RP)/$(main) $(OP)/$(main).o $(OP)/SIAM.o $(OP)/Grid.o $(OP)/log.o $(OP)/Params.o $(OP)/routines.o $(OP)/dinterpl.o $(OP)/SIAM_GPU.o $(LIBS) $(CUDALIBS)

# main program
$(OP)/$(main).o : $(SP)/$(main).cpp $(SP)/SIAM.h $(SP)/Grid.h
	$(mpiCC) $(FLAGS) -c -o $@ $(SP)/$(main).cpp

# SIAM
$(OP)/SIAM.o : $(SP)/SIAM.cpp $(SP)/SIAM.h $(SP)/Grid.h $(SP)/Params.h $(SP)/routines.h $(SP)/log.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/SIAM.cpp

# cuSIAM (GPU)
$(OP)/SIAM_GPU.o : $(SP)/SIAM.cu $(SP)/SIAM.h
	nvcc $(CUDAFLAGS) -c -o $@ $(SP)/SIAM.cu

# Result
$(OP)/Grid.o : $(SP)/Grid.cpp $(SP)/Grid.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Grid.cpp
	
# Logging
$(OP)/log.o : $(SP)/log.cpp $(SP)/log.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/log.cpp

# Input class used for reading files with parameters
$(OP)/Params.o : $(SP)/Params.cpp $(SP)/Params.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Params.cpp

# Interpolation
$(OP)/dinterpl.o : $(SP)/dinterpl.cpp $(SP)/dinterpl.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/dinterpl.cpp

# contains some constants and useful numerical routines
$(OP)/routines.o : $(SP)/routines.cpp $(SP)/routines.h 
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/routines.cpp

# clean all object and exec files
clean :
	rm -f $(RP)/$(main) $(OP)/*.o
