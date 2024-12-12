CXX=g++
CXXFLAGS = -std=c++11 -O3 -Wall -g
LIBS= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lm -lboost_program_options

CXX_SOURCES= matrix_utils.cpp read_parameters.cpp file.cpp regem.cpp center.cpp
CXX_OBJECTS= matrix_utils.o read_parameters.o file.o regem.o center.o
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
