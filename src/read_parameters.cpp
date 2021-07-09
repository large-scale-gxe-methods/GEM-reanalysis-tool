/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/


#include "declars.h"
#define VERSION "1.0"

void print_help();


// Function to process command line arguments
void CommandLine::processCommandLine(int argc, char* argv[]) {


    cout << "\n*********************************************************\n";
    cout << "Welcome to ReGEM v" << VERSION << "\n";
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
        ("out", po::value<std::string>()->default_value("gem_reanalysis.out"), "")
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
    if (out.count("output-style")) {
        outStyle = out["output-style"].as<string>();
        if (outStyle.compare("minimum") != 0 && outStyle.compare("meta") != 0  && outStyle.compare("full") != 0 ) {
            cerr << "\nERROR: --output-style should be minimum, meta or full.\n\n";
            exit(1);
        }
    }

    // Exposure
    if (out.count("exposure-names")) {
        exp = out["exposure-names"].as<std::vector<std::string>>();
        if (exp.size() == 0) {
            cerr << "\nERROR: No exposure (--exposure-names) is specified.\n\n";
        }
        std::unordered_map<std::string, int> expHM;
        for (size_t i = 0; i < exp.size(); i++) {
            expHM[exp[i]] += 1;
            if (expHM[exp[i]] > 1) {
                cerr << "\nERROR: Exposure " + exp[i] + " is specified more than once.\n\n";
                exit(1);
            }
        }
        numExpSelCol = exp.size();
    }
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

    
    // Print parameter info
    cout << "The Input File is: " << inFile << "\n\n";

    if (numIntSelCol == 0) {
        cout << "No Interaction Covariates Selected" << "\n";
    }
    else {
        cout << "The Total Number of Selected Interaction Covariates is: " << numIntSelCol << "\n";
        cout << "The Selected Interaction Covariates are:  ";
        for (int i = 0; i < numIntSelCol; i++) {
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
            cout << exp[i] << "   ";
        }
        cout << "\n\n";
    }

    cout << "Output File: " << outFile << "\n";
    cout << "*********************************************************\n";

}






void print_help() {

 
    cout << "General Options: " << endl
        << "   --help \t\t Prints available options and exits." << endl
        << "   --version \t\t Prints the version of GEM and exits." << endl;
    cout << endl << endl;



    cout << "Input File Options: " << endl
        << "   --pheno-file \t Path to the phenotype file." << endl
        << "   --bgen \t\t Path to the BGEN file." << endl
        << "   --sample \t\t Path to the sample file. Required when the BGEN file does not contain sample identifiers." << endl
        << "   --pfile \t\t Path and prefix to the .pgen, .pvar, and .psam files." << endl
        << "   --pgen \t\t Path to the pgen file." << endl
        << "   --pvar \t\t Path to the pvar file." << endl
        << "   --psam \t\t Path to the psam file." << endl
        << "   --bfile \t\t Path and prefix to the .bed, .bim and .fam files." << endl
        << "   --bed \t\t Path to the bed file." << endl
        << "   --bim \t\t Path to the bim file." << endl
        << "   --fam \t\t Path to the fam file." << endl
        << "   --out \t\t Full path and extension to where GEM output results. \n \t\t\t    Default: gem.out" << endl;
    cout << endl << endl;


    cout << "Phenotype File Options: " << endl
        << "   --sampleid-name \t Column name in the phenotype file that contains sample identifiers." << endl
        << "   --pheno-name \t Column name in the phenotype file that contains the phenotype of interest." << endl
        << "   --exposure-names \t One or more column names in the phenotype file naming the exposure(s) to be included in interaction tests." << endl
        << "   --int-covar-names \t Any column names in the phenotype file naming the covariate(s) for which interactions should\n \t\t\t   be included for adjustment (mutually exclusive with --exposure-names)." << endl
        << "   --covar-names \t Any column names in the phenotype file naming the covariates for which only main effects should\n \t\t\t   be included for adjustment (mutually exclusive with both --exposure-names and --int-covar-names)." << endl
        << "   --pheno-type \t 0 indicates a continuous phenotype and 1 indicates a binary phenotype." << endl
        << "   --robust \t\t 0 for model-based standard errors and 1 for robust standard errors. \n \t\t\t    Default: 0" << endl
        << "   --tol \t\t Convergence tolerance for logistic regression. \n \t\t\t    Default: 0.0000001" << endl
        << "   --delim \t\t Delimiter separating values in the phenotype file. Tab delimiter should be represented as \\t and space delimiter as \\0. \n \t\t\t    Default: , (comma-separated)" << endl
        << "   --missing-value \t Indicates how missing values in the phenotype file are stored. \n \t\t\t    Default: NA" << endl
        << "   --center \t\t 0 for no centering to be done and 1 to center ALL exposures and covariates. \n \t\t\t    Default: 1" << endl
        << "   --scale \t\t 0 for no scaling to be done and 1 to scale ALL exposures and covariates by the standard deviation. \n \t\t\t    Default: 0" << endl;
    cout << endl << endl;


    cout << "Filtering Options: " << endl
        << "   --maf \t\t Threshold to filter variants based on the minor allele frequency.\n \t\t\t    Default: 0.001" << endl
        << "   --miss-geno-cutoff \t Threshold to filter variants based on the missing genotype rate.\n \t\t\t    Default: 0.05" << endl
        << "   --include-snp-file \t Path to file containing a subset of variants in the specified genotype file to be used for analysis. The first\n \t\t\t   line in this file is the header that specifies which variant identifier in the genotype file is used for ID\n \t\t\t   matching. This must be 'snpid' or 'rsid' (BGEN only). There should be one variant identifier per line after the header.\n \t\t\t   Variants not listed in this file will be excluded from analysis." << endl;
    cout << endl << endl;



    cout << "Performance Options:" << endl
        << "   --threads \t\t Set number of compute threads \n \t\t\t    Default: ceiling(detected threads / 2)" << endl
        << "   --stream-snps \t Number of SNPs to analyze in a batch. Memory consumption will increase for larger values of stream-snps.\n \t\t\t    Default: 1" << endl;
    cout << endl << endl;

}