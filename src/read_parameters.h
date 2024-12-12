#ifndef READPARAMETERS2_H
#define READPARAMETERS2_H

#include <vector>
#include <string>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>

class CommandLine {

public:
	// Phenotype file
	std::string inFile;
	std::string delim;
	char pheno_delim;
	
	// Out file
	std::string outFile;
	std::string outStyle; 

	std::vector<std::string> interactions;

	size_t nExp  = 0;
	size_t nExp1 = 0;
	size_t nIntCov = 0;
	size_t nInt1 = 0;

	// Ouput style
	int  printStart = 0;
	int  printEnd   = 0;

	// Center conversion options
	int centerIn = 0;
	int centerOut = 0;
	std::unordered_map<std::string, double> meanValues;

	// Performance options
	int threads;
	void processCommandLine(int argc, char* argv[]);
};

void print_welcome();

void print_help();

#endif
