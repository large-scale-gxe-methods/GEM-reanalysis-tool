/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/

#include "declars.h"

#define VERSION "1.0"

void print_help();

// Function to process command line arguments
void CommandLine::processCommandLine(int argc, char* argv[]) {


    cout << "\n*********************************************************\n";
    cout << "Welcome to REGEM v" << VERSION << "\n";
    cout << "(C) 2021 Duy Pham and Han Chen \n";
    cout << "GNU General Public License v3\n";
    cout << "*********************************************************\n";


    // GEM parameters. Details are printed from the print_help() function below.

    // General commands
    po::options_description general("General options");
    general.add_options()
        ("help, h", "")
        ("version", "");

    // Input/Output file options
    po::options_description files("Input/Output file options");
    files.add_options()
        ("input-file", po::value<std::string>(), "")
        ("out", po::value<std::string>()->default_value("regem.out"), "")
        ("output-style", po::value<std::string>()->default_value("minimum"), "");

    // Phenotype file
    po::options_description resultsfile("Results file options");
    resultsfile.add_options()
        ("exposure-names", po::value<std::vector<std::string>>()->multitoken(), "")
        ("int-covar-names", po::value<std::vector<std::string>>()->multitoken(), "");

    // Combine all options together
    po::options_description all("Options");
    all.add(general).add(files).add(resultsfile);

    po::variables_map out;


    // QC
    try {
        po::store(po::command_line_parser(argc, argv)
            .options(all)
            .style(po::command_line_style::unix_style
                | po::command_line_style::allow_long_disguise)
            .run(),
            out);

    }
    catch (po::error const& e) {
        std::cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    po::notify(out);



    // General
    if (out.count("help")) {
        print_help();
        exit(1);

    }
    if (out.count("version")) {
        cout << "\nTest version: " << VERSION << "\n" << endl;
        exit(1);

    }


    // Input/Output Files
    if (out.count("input-file")) {
        inFile = out["input-file"].as<std::string>();

    }
    else {
        cerr << "\nERROR: Results file (--input-file) is needed. \n\n";
        exit(1);

    }


    // Interaction covariates
    if (out.count("int-covar-names")) {
        icov = out["int-covar-names"].as<std::vector<std::string>>();

        std::unordered_map<std::string, int> icovHM;
        for (size_t i = 0; i <  icov.size(); i++) {
            icovHM[icov[i]] += 1;
            if (icovHM[icov[i]] > 1) {
                cerr << "\nERROR: Interaction covariate " + icov[i] + " is specified more than once.\n\n";
                exit(1);
            }
        }
        numIntSelCol = icov.size();
    }

    // Exposures
    if (out.count("exposure-names")) {
        interactions = out["exposure-names"].as<std::vector<std::string>>();
        if (interactions.size() == 0) {
            cerr << "\nERROR: No exposure (--exposure-names) is specified.\n\n";
        }
        std::unordered_map<std::string, int> expHM;
        for (size_t i = 0; i < interactions.size(); i++) {
            expHM[interactions[i]] += 1;
            if (expHM[interactions[i]] > 1) {
                cerr << "\nERROR: Exposure " + interactions[i] + " is specified more than once.\n\n";
                exit(1);
            }
        }

        numExpSelCol = interactions.size() - numIntSelCol;
        if (numIntSelCol > 0) {
		    for (int i = 0; i < numIntSelCol; i++)
			    interactions.push_back(icov[i]);
        }
    }


    // Output file
    if (out.count("out")) {
        outFile = out["out"].as<string>();

        std::ofstream results(outFile, std::ofstream::binary);
        if (!results) {
            cerr << "\nERROR: Output file could not be opened.\n\n";
            exit(1);
        }

        if (results.fail()) {
            cerr << "\nERROR: Output file could not be opened.\n\n";
            exit(1);
        }

        results << "test" << endl;
        if (results.fail()) {
            cerr << "\nERROR: Cannot write to output file.\n\n";
            results.close();
            boost::filesystem::remove(outFile.c_str());
            exit(1);
        }
        results.close();

        boost::filesystem::remove(outFile.c_str());
    }



    // Output style
    if (out.count("output-style")) {

        outStyle = out["output-style"].as<string>();
        if (outStyle.compare("minimum") != 0 && outStyle.compare("meta") != 0  && outStyle.compare("full") != 0 ) {
            cerr << "\nERROR: --output-style should be minimum, meta or full.\n\n";
            exit(1);
        }

        if (outStyle.compare("meta") == 0) {
            printStart = 0; 
            printEnd   = numIntSelCol + numExpSelCol + 1; 
            printMeta  = true;
        } else if (outStyle.compare("full") == 0) {
            printStart = 0; 
            printEnd   = numIntSelCol + numExpSelCol + 1; 
            printFull  = true;
        } else {
            printStart = 1;
            printEnd   = numExpSelCol + 1;
        }
    }



    // Print parameter info
    cout << "The Input File is: " << inFile << "\n\n";

    if (numIntSelCol == 0) {
        cout << "No Interaction Covariates Selected" << "\n";
    }
    else {
        cout << "The Total Number of Selected Interaction Covariates is: " << numIntSelCol << "\n";
        cout << "The Selected Interaction Covariates are:  ";
        for (int i = numExpSelCol; i < numIntSelCol+1; i++) {
            cout << icov[i] << "   ";
        }
        cout << "\n";
    }

    if (numExpSelCol == 0) {
        cout << "No Exposures Selected" << "\n";
    }
    else {
        cout << "The Total Number of Exposures is: " << numExpSelCol << '\n';

        cout << "The Selected Exposures are:  ";
        for (int i = 0; i < numExpSelCol; i++) {
            cout << interactions[i] << "   ";
        }
        cout << "\n\n";
    }

    cout << "Output File: " << outFile << "\n";
    cout << "*********************************************************\n";

}






void print_help() {

 
    cout << "General Options: " << endl
        << "   --help \t\t Prints available options and exits." << endl
        << "   --version \t\t Prints the version of REGEM and exits." << endl;
    cout << endl << endl;



    cout << "Input File Options: " << endl
        << "   --results-file \t Path to the GEM results file." << endl
        << "   --out \t\t Full path and extension to where REGEM output results. \n \t\t\t    Default: regem.out" << endl;
    cout << endl << endl;


    cout << "Phenotype File Options: " << endl
        << "   --exposure-names \t One or more column names in the phenotype file naming the exposure(s) to be included in interaction tests." << endl
        << "   --int-covar-names \t Any column names in the phenotype file naming the covariate(s) for which interactions should\n \t\t\t   be included for adjustment (mutually exclusive with --exposure-names)." << endl;
    cout << endl << endl;


    cout << "Filtering Options: " << endl
        << "   --maf \t\t Threshold to filter variants based on the minor allele frequency.\n \t\t\t    Default: 0.001" << endl
        << "   --miss-geno-cutoff \t Threshold to filter variants based on the missing genotype rate.\n \t\t\t    Default: 0.05" << endl;
    cout << endl << endl;

    cout << endl << endl;

}