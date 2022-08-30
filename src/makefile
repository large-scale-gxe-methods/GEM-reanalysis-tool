CXX=g++
CXXFLAGS = -std=c++11 -O3 -Wall -g

# Uncomment line below for LAPACK/BLAS only and comment other lines beginning with LIBS.
#LIBS= -llapack -lblas -ldl -pthread -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem

LIBS= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lm -lboost_program_options

CXX_SOURCES= matrix_utils.cpp read_parameters.cpp file.cpp regem.cpp
CXX_OBJECTS= matrix_utils.o read_parameters.o file.o regem.o
OBJS= $(CXX_OBJECTS)

#targets
all: REGEM

REGEM: $(CXX_OBJECTS) $(C_OBJECTS)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS)
