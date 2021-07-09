CC=gcc  # Do not change this for libdeflate package.
CXX=g++
CXXFLAGS = -std=c++11 -O2 -g

# Uncomment line below for LAPACK/BLAS only and comment other lines beginning with LIBS.
#LIBS= -llapack -lblas -ldl -pthread -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem

LIBS= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem

CXX_SOURCES= matrix_utils.cpp read_parameters.cpp regem_main.cpp
CXX_OBJECTS= matrix_utils.o read_parameters.o regem_main.o
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
