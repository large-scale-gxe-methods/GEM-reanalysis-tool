#ifndef REGEM_H
#define REGEM_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <numeric>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include "matrix_utils.h"
#include "read_parameters.h"
#include "file.h"
#include "center.h"

#define VERSION "1.0"

using std::cout;
using std::endl;
using std::cerr;

void regem(CommandLine* cmd);

void printTimeCompleted(double wall0, double wall1, double cpu0, double cpu1);

#endif
