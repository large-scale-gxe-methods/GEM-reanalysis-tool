#ifndef READPARAMETERS2_H
#define READPARAMETERS2_H

#include <string>
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

	std::vector<std::string> icov;	
	std::vector<std::string> interactions;

	int numExpSelCol = 0;
	int numIntSelCol = 0;

	// Ouput style
	bool printMeta  = false;
	bool printFull  = false;
	int  printStart;
	int  printEnd;

	// Performance options
	int threads;
	void processCommandLine(int argc, char* argv[]);
};


#endif
