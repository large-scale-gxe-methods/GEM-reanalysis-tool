/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/

#include "regem.h"


namespace po = boost::program_options;

// Function to process command line arguments
void CommandLine::processCommandLine(int argc, char* argv[]) 
{
    print_welcome();

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
        ("output-style", po::value<std::string>()->default_value("full"), "");

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
        std::cout << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    po::notify(out);


    // General
    if (out.count("help")) {
        print_help();
        exit(1);

    }
    if (out.count("version")) {
        cout << "\nREGEM version: " << VERSION << "\n" << endl;
        exit(1);

    }


    // Input/Output Files
    if (out.count("input-file")) {
        inFile = out["input-file"].as<std::string>();

    }
    else {
        cout << "\nERROR: Results file (--input-file) is needed. \n\n";
        exit(1);

    }


    // Exposures
    if (out.count("exposure-names")) {
        interactions = out["exposure-names"].as<std::vector<std::string>>();

        std::set<std::string> s(interactions.begin(), interactions.end());
        if (s.size() != interactions.size()) {
            cout << "\nERROR: There are duplicate exposure names (--exposure-names).\n\n";
            exit(1);
        }

        nExp  = interactions.size();
        nExp1 = nExp + 1;
    } else {
        cout << "\nERROR: No exposure (--exposure-names) specified.\n\n";
        exit(1);
    }

    // Interaction covariates
    if (out.count("int-covar-names")) {
        std::vector<std::string> icov = out["int-covar-names"].as<std::vector<std::string>>();

        std::set<std::string> s(icov.begin(), icov.end());
        if (s.size() != icov.size()) {
            cout << "\nERROR: There are duplicate interaction covariates names (--int-cov-names).\n\n";
            exit(1);
        }

        nIntCov = icov.size();
        for (size_t i = 0; i < nIntCov; i++) {
            if (std::find(interactions.begin(), interactions.end(), icov[i]) != interactions.end()) {
                cout << "\nERROR: Interaction covariate " + icov[i] + " is already specified in --exposure-names.\n\n";
            }
            interactions.push_back(icov[i]);
        }
    }
    nInt1 = interactions.size() + 1;


    // Output file
    if (out.count("out")) {
        outFile = out["out"].as<std::string>();

        std::ofstream results(outFile, std::ofstream::binary);
        if (!results) {
            cout << "\nERROR: Output file could not be opened.\n\n";
            exit(1);
        }

        if (results.fail()) {
            cout << "\nERROR: Output file could not be opened.\n\n";
            exit(1);
        }

        results << "test" << endl;
        if (results.fail()) {
            cerr << "\nERROR: Cannot write to output file.\n\n";
            results.close();
            
            if (std::remove(outFile.c_str()) != 0) {
                cerr << "\nERROR: Cannot delete output file.\n\n";
            }
            exit(1);
        }
        results.close();

        if (std::remove(outFile.c_str()) != 0) {
            cerr << "\nERROR: Cannot delete output file.\n\n";
            exit(1);
        }
    }


    // Output style
    if (out.count("output-style")) {
        outStyle = out["output-style"].as<std::string>();
        if (outStyle.compare("minimum") != 0 && outStyle.compare("meta") != 0  && outStyle.compare("full") != 0 ) {
            cout << "\nERROR: --output-style should be minimum, meta or full.\n\n";
            exit(1);
        }

        if (outStyle.compare("meta") == 0 || outStyle.compare("full") == 0) {
            printStart = 0; 
            printEnd   = nInt1; 
        } else {
            printStart = 1;
            printEnd   = nExp1;
        }
    }


    // Print parameter info
    cout << "The input file is [" << inFile << "].\n\n";
    if (nIntCov == 0) {
        cout << "No interaction covariates selected." << "\n";
    }
    else {
        cout << "The selected interaction covariate(s) are: ";
        for (size_t i = nExp; i < (nInt1 - 1); i++) {
            cout << interactions[i] << " ";
        }
        cout << "\n";
    }
    cout << "The selected exposure(s) are: ";
    for (size_t i = 0; i < nExp; i++) {
        cout << interactions[i] << "  ";
    }
    cout << "\n";
    cout << "*********************************************************\n";
}


void print_welcome() {
    cout << "\n*********************************************************\n";
    cout << "Welcome to REGEM v" << VERSION << "\n";
    cout << "(C) 2021-2022 Duy Pham and Han Chen \n";
    cout << "GNU General Public License v3\n";
    cout << "*********************************************************\n";
}

void print_help() {
    cout << "General Options: " << endl
        << "   --help \t\t Prints available options and exits." << endl
        << "   --version \t\t Prints the version of REGEM and exits." << endl;
    cout << endl << endl;



    cout << "Input File Options: " << endl
        << "   --input-file \t Path to the GEM results file." << endl
        << "   --out \t\t Full path and extension to where REGEM output results. \n \t\t\t    Default: regem.out" << endl
        << "   --output-style \t Modifies the output of REGEM. Must be one of the following: \n\t\t\t    minimum: Output the summary statistics for only the GxE and marginal G terms. \n \t\t\t    meta: 'minimum' output plus additional fields for the main G and any GxCovariate terms \n \t\t\t\t  For a robust analysis, additional columns for the model-based summary statistics will be included.  \n \t\t\t    full: 'meta' output plus additional fields needed for re-analyses of a subset of interactions \n \t\t\t    Default: full" << endl; 
    cout << endl << endl;


    cout << "Phenotype File Options: " << endl
        << "   --exposure-names \t One or more column names in the phenotype file naming the exposure(s) to be included in interaction tests." << endl
        << "   --int-covar-names \t Any column names in the phenotype file naming the covariate(s) for which interactions should\n \t\t\t   be included for adjustment (mutually exclusive with --exposure-names)." << endl;
    cout << endl << endl;
    cout << endl << endl;

}
