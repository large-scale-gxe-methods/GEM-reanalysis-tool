
#include "declars.h"

static inline double Square (double x) { return x*x; }
void printExecutionTime(std::chrono::high_resolution_clock::time_point start_time, std::chrono::high_resolution_clock::time_point end_time);
void printOutputHeader(int n_varInfo, int printStart, int printEnd, int robust, bool strata, vector<string> catE_headers, vector<string> covNames, string outStyle, string output);
std::vector<std::string> cartesian_vec( vector<vector<string> >& v );
std::vector<std::string> cartesian_vec_sep( vector<vector<string> >& v);

int main(int argc, char* argv[]) 
{
	CommandLine cmd;
	cmd.processCommandLine(argc, argv);
	auto start_time = std::chrono::high_resolution_clock::now();
	string inFile  = cmd.inFile;
	string outFile = cmd.outFile;
	vector<string> exp  = cmd.exp;
	vector<string> icov = cmd.icov;
	int expSq  = cmd.numExpSelCol;
	int expSq1 = expSq+1;
	int intSq  = cmd.numIntSelCol;
	int Sq1    = intSq + expSq1;
	int robust = 0;
	double sigma2 = 0.0;
	bool strata = false;
	boost::math::chi_squared chisq_dist_M(1);
	boost::math::chi_squared chisq_dist_Int(expSq);
	boost::math::chi_squared chisq_dist_Joint(expSq1);

    int printStart = 1; 
	int printEnd   = expSq1;
    bool printFull = false;
	bool printMeta = false;
    if (cmd.outStyle.compare("meta") == 0) {
        printStart = 0; printEnd = Sq1; printMeta = true;
    } else if (cmd.outStyle.compare("full") == 0) {
        printStart = 0; printEnd = Sq1; printFull = true;
    }


	if (intSq > 0 ) {
		for (int i = 0; i < intSq; i++)
			exp.push_back(icov[i]);
	}


	// Read input file
	std::ifstream res;
    res.open(inFile);
    if (!res.is_open()) {
        cerr << "\nERROR: Cannot open results file. \n\n" << endl;
        exit(1);
    }

	// Get dispersion
	std::string line;
    getline(res, line);
	if (line.rfind("#dispersion", 0) != 0) {
		cerr << "\nERROR: The first line in the input file is not '#dispersion'.\n\n";
		exit(1); 
	} else {
		std::istringstream iss(line);
		string value;
		vector <string> values;
		while (getline(iss, value, ' ')) values.push_back(value);
		sscanf(values[1].c_str(), "%lf", &sigma2);
	}


	// Read the column header
	getline(res, line);
    std::istringstream issHead(line);
    string header;
    int header_i = 0;
	std::unordered_map<std::string, int> colNames;
	vector<string> strata_names;
	int n_strata = 0;
    while (getline(issHead, header, '\t')) {
        header.erase(std::remove(header.begin(), header.end(), '\r'), header.end());
		if (colNames.find(header) != colNames.end()) {
            cerr << "\nERROR: There are duplicate header names (" << header << ") in the results file.\n\n";
            exit(1);
        }
		colNames[header] = header_i;

        if (header.rfind("Beta_G-", 0) == 0) {
            string tmp_name = header;
            tmp_name.erase(0, 7);
            if (std::find(exp.begin(), exp.end(), tmp_name) == exp.end()) {
                exp.push_back(tmp_name);
            }
        }

		if (header.rfind("N_", 0) == 0) {
			if (header != "N_Samples") {
				strata_names.push_back(header);
				n_strata++;
			}
		}

		if ((robust !=1) && (header.rfind("robust_", 0) == 0)) {
			robust = 1;
		}
        header_i++;
    }
	int n_varInfo = colNames["AF"];


	std::vector<int> rgs_len;
	std::vector<int> rgs_index;
	std::vector<string> catE_headers;
	std::vector<vector<string>> rg_strata;
	size_t n_cart = 0;
	if (strata_names.size() > 0) {
		std::string line2;
		std::vector<std::string> vec;
		std::stringstream ss(strata_names[0]);
		int vecLen = 0;
		while(std::getline(ss, line2, '_')) {
			vec.push_back(line2);
			vecLen++;
		}

		size_t vs = vec.size();
		int k = 1;
		string tmp = "";
		vector<string> catE;
		int n_catE = 0;
		while (k < vs) {
			tmp = (tmp.length() == 0) ? tmp + vec[k] : tmp + "_" + vec[k];
			if (std::find(exp.begin(), exp.end(), tmp) != exp.end()) {
				catE.push_back(tmp);
				k++;
				vs--;
				tmp = "";
				n_catE++;
			}
		}


		int n_rgs = 0;
		vector<int> rg_si;
		std::vector<string> rg_catE;
		for (int i = 0; i < (expSq + intSq); i++) {
			auto it = std::find(catE.begin(), catE.end(), exp[i]);
			if (it != catE.end()) {
				int index = it - catE.begin();
				rg_si.push_back(vecLen - (n_catE - index));
				rg_catE.push_back(exp[i]);
				n_rgs++;
			}
		}

		if (n_rgs > 0 ) {
			strata = true;
			std::vector<std::string> vec;
			rg_strata.resize(n_rgs);
			for (size_t i = 0; i < strata_names.size(); i++) {
				std::string line2;
				std::stringstream ss(strata_names[i]);
				int vecLen = 0;
				while(std::getline(ss, line2, '_')) {
					vec.push_back(line2);
					vecLen++;
				}

				for (int j = 0; j < n_rgs; j++) {
					if (std::find(rg_strata[j].begin(), rg_strata[j].end(), vec[rg_si[j]]) == rg_strata[j].end()) {
						rg_strata[j].push_back(vec[rg_si[j]]);
					}
				}
				vec.clear();
			}

			std::vector<std::string> rgs_vec = cartesian_vec(rg_strata);
			n_cart = rgs_vec.size();

			std::vector<std::string> vec2;
			for(size_t i = 0; i < rgs_vec.size(); i++) {
				int cnt = 0;
				for (size_t j = 0; j < strata_names.size(); j++){
					std::string line2;
					std::stringstream ss(strata_names[j]);
					while(std::getline(ss, line2, '_')) {
						vec2.push_back(line2);
					}

					string tmp = "";
					for (int k = 0; k < n_rgs; k++) {
						tmp += vec2[rg_si[k]];
					}
					if (tmp == rgs_vec[i]) {
						rgs_index.push_back(j*2);
						cnt++;
					} 
					vec2.clear();
				}
				rgs_len.push_back(cnt);
			}

			string catE = "";
			for (int i = 0; i < n_rgs; i++)
				catE += rg_catE[i] + "_";
			std::vector<std::string> rgs_cart_sep = cartesian_vec_sep(rg_strata);
			for (size_t i = 0; i < rgs_cart_sep.size(); i++) {
				catE_headers.push_back(catE + rgs_cart_sep[i]);
			}

		}	
	} 
	

	// Get column header index for marginal, G, GxE, and GxC
	int beta_marginal   = colNames["Beta_Marginal"];
	int M_se_marginal   = colNames["SE_Beta_Marginal"];
	int R_se_marginal   = colNames["robust_SE_Beta_Marginal"];
	int M_pval_marginal = colNames["P_Value_Marginal"];
	int R_pval_marginal = colNames["robust_P_Value_Marginal"];

	for (size_t i = 0; i < exp.size(); i++) {
		exp[i] = "G-" + exp[i];
	}
	exp.insert(exp.begin(), "G");

	int dim = exp.size();
	vector<int> beta_idx(dim);
	vector<int> mb_idx(dim*dim);	
	vector<int> rb_idx(dim*dim);
	for (int i = 0; i < dim; i++) {
		// Beta columns
		if (colNames.find("Beta_" + exp[i]) == colNames.end()) {
			cerr << "\nERROR: Column header Beta_" + exp[i] + " does not exist in input file.\n\n";
			exit(1);
		} else {
			beta_idx[i] = colNames["Beta_" + exp[i]];
		}

		// Model-based SE columns
		if (colNames.find("SE_Beta_" + exp[i]) == colNames.end()) {
			cerr << "\nERROR: Column header SE_Beta_" + exp[i] + " does not exist in the input file.\n\n";
			exit(1);
		} else {
			mb_idx[(i * dim) + i] = colNames["SE_Beta_" + exp[i]];	
		}

		// Robust SE columns
		if (robust == 1) {
			if (colNames.find("robust_SE_Beta_" + exp[i]) == colNames.end()) {
				cerr << "\nERROR: Column header robust_SE_Beta_" + exp[i] + " does not exist in the input file.\n\n";
				exit(1);
			} else {
				rb_idx[(i * dim) + i] = colNames["robust_SE_Beta_" + exp[i]];	
			}
		}

		// Model-based Covariance columns
		for (int j = (i+1); j < dim; j++) {
			if (colNames.find("Cov_Beta_" + exp[i] + "_" + exp[j]) != colNames.end()) {
				int tmp_idx = colNames["Cov_Beta_" + exp[i] + "_" + exp[j]];
				mb_idx[(i * dim) + j] = tmp_idx;
				mb_idx[(j * dim) + i] = tmp_idx;
			} 
			else if (colNames.find("Cov_Beta_" + exp[j] + "_" + exp[i]) != colNames.end()) {
 				int tmp_idx = colNames["Cov_Beta_" + exp[j] + "_" + exp[i]];
				mb_idx[(i * dim) + j] = tmp_idx;
				mb_idx[(j * dim) + i] = tmp_idx;
			} else {
				cerr << "\nERROR: Column header Cov_Beta_" + exp[i] + "_" + exp[j] + " does not exist in the input file.\n\n";
				exit(1);
			}
		}

		// Robust Covariance columns
		if (robust == 1) {
			for (int j = (i+1); j < dim; j++) {
				if (colNames.find("robust_Cov_Beta_" + exp[i] + "_" + exp[j]) != colNames.end()) {
					int tmp_idx = colNames["robust_Cov_Beta_" + exp[i] + "_" + exp[j]];
					rb_idx[(i * dim) + j] = tmp_idx;
					rb_idx[(j * dim) + i] = tmp_idx;
				}
				else if (colNames.find("robust_Cov_Beta_" + exp[j] + "_" + exp[i]) != colNames.end()) {
					int tmp_idx = colNames["robust_Cov_Beta_" + exp[j] + "_" + exp[i]];
					rb_idx[(i * dim) + j] = tmp_idx;
					rb_idx[(j * dim) + i] = tmp_idx;					
				}
				else {
					cerr << "\nERROR: Column header Cov_Beta_" + exp[i] + "_" + exp[j] + " does not exist in the input file.\n\n";
					exit(1);
				}
			}	
		}
	}



	// Check if there are results in the input file
    int nvars = 0;
    while (getline(res, line)) nvars++;
	if (nvars == 0) {
		cerr << "\nERROR: No variants in the input file.\n\n.";
		exit(1);
	} else {
		res.clear(); 
		res.seekg(0, res.beg);
    	getline(res, line);
		getline(res, line);	
	}

	// Create the output file and write column header names
	std::ofstream results(outFile, std::ios_base::app);
	std::ostringstream oss;
	printOutputHeader(n_varInfo, printStart, printEnd, robust, strata, catE_headers, exp, cmd.outStyle, outFile);



	// REGEM begins
	vector<double> beta_vec(dim);
	vector<double> mb_vec(dim*dim);
	vector<double> rb_vec(dim*dim);
	vector<double> strata_N(n_cart);
	vector<double> strata_AF(n_cart);
	double* Beta = &beta_vec[0];
	double* MV   = &mb_vec[0];
	double* RV   = &rb_vec[0];

	double* StempE   = new double[expSq];
	double* StempGE  = new double[expSq1];
	double* U        = new double[dim];
	double* Beta_dot = new double[Sq1];
	double* V        = new double[dim * dim];
	double* Vi_dot    = new double[Sq1*Sq1];
	double* VE_dot   = new double[expSq*expSq];
	double* O        = new double[dim*dim];
	double* O_dot    = new double[Sq1*Sq1];
	double* RV_dot   = new double[Sq1*Sq1];
	double* A_i      = new double[expSq1 * expSq1];
	double M_pvalM, M_pvalInt, M_pvalJoint, R_pvalM, R_pvalInt, R_pvalJoint;

	int k = 0;
	for (int p = 0; p < nvars; p++) {
		getline(res, line);
		std::istringstream iss(line);
		string value;
		vector <string> values;
		while (getline(iss, value, '\t')) values.push_back(value);

		for (int i = 0; i < dim; i++) {
			sscanf(values[beta_idx[i]].c_str(), "%lf", &beta_vec[i]);
			for (int j = 0; j < dim; j++) {
	 			sscanf(values[mb_idx[(i * dim) + j]].c_str(), "%lf", &mb_vec[(i * dim) + j]);
	 		}
			mb_vec[(i * dim) + i] *= mb_vec[(i * dim) + i];
		}
		if (robust == 1) {
			for (int i = 0; i < dim; i++) {
				for (int j = 0; j < dim; j++) {
					sscanf(values[rb_idx[(i * dim) + j]].c_str(), "%lf", &rb_vec[(i * dim) + j]);
				}
				rb_vec[(i * dim) + i] *= rb_vec[(i * dim) + i];
			}
		}

		// Compute V from model-based variance-covariance matrix
		subMatrix(MV, V, dim, dim, dim, dim, 0);
		matInv(V, dim);
		for (int i = 0; i < dim*dim; i++)
			V[i] *= sigma2;

		subMatrix(V, Vi_dot, Sq1, Sq1, dim, Sq1, 0);
		matInv(Vi_dot, Sq1);

		// Compute U
 		matvecprod(V, Beta, U, dim, dim);

		// Compute new betas
		matvecprod(Vi_dot, U, Beta_dot, Sq1, Sq1);

	
		if (robust == 1) {
			//Compute Omega
			matNmatNprod(V, RV, O, dim, dim, dim);
	 		matNmatNprod(O, V,  O, dim, dim, dim);
			subMatrix(O, O_dot, Sq1, Sq1, dim, Sq1, 0);

			// Compute robust variance-covariance
            matNmatNprod(O_dot, Vi_dot,  RV_dot, Sq1, Sq1, Sq1);
			matNmatNprod(Vi_dot, RV_dot, RV_dot, Sq1, Sq1, Sq1);

			// Robust Stat Int
			subMatrix(RV_dot, VE_dot, expSq, expSq, Sq1, expSq, Sq1 + 1);
			matInv(VE_dot, expSq);
			matvecSprod(VE_dot, Beta_dot, StempE, expSq, expSq, 1);

			double statInt = 0.0;
			for (int j = 1; j < expSq1; j++) 
				statInt += Beta_dot[j] * StempE[j-1];

			R_pvalInt = (isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));

			// Robust Stat Joint
			subMatrix(RV_dot, A_i, expSq1, expSq1, Sq1, expSq1, 0);
			matInv(A_i, expSq1);
			matvecprod(A_i, Beta_dot, StempGE, expSq1, expSq1);

			double statJoint = 0.0;
			for (int k = 0; k < expSq1; k++)
				statJoint += Beta_dot[k] * StempGE[k];
			
			R_pvalJoint = (isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));
		}


		// Get model-based variance-covariance 
		for (int i = 0; i < Sq1 * Sq1; i++) {
			Vi_dot[i] *= sigma2;
		}

		// Model-based Stat Int
		subMatrix(Vi_dot, VE_dot, expSq, expSq, Sq1, expSq, Sq1 + 1);
		matInv(VE_dot, expSq);
		matvecSprod(VE_dot, Beta_dot, StempE, expSq, expSq, 1);

		double statInt = 0.0;
		for (int j = 1; j < expSq1; j++) 
			statInt += Beta_dot[j] * StempE[j-1];

		M_pvalInt = (isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));

		// Model-based Stat Joint
		subMatrix(Vi_dot, A_i, expSq1, expSq1, Sq1, expSq1, 0);
		matInv(A_i, expSq1);
		matvecprod(A_i, Beta_dot, StempGE, expSq1, expSq1);

		double statJoint = 0.0;
		for (int k = 0; k < expSq1; k++)
			statJoint += Beta_dot[k] * StempGE[k];
			
		M_pvalJoint = (isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));

		if (strata) {
			size_t i = 0;
			size_t k = 0;

			vector<double> tv(n_strata*2);
			for (int j = 0; j < (n_strata*2); j++) {
				tv[j] = std::stod(values[(n_varInfo+1) + j]);
			}

            while (i < rgs_index.size()) {
                for (int j = 0; j < rgs_len[k]; j++) {
                    strata_N[k]  += tv[rgs_index[i]];
                    strata_AF[k] += (tv[rgs_index[i]] * tv[rgs_index[i]+1]);
                    i++;
                }
                strata_AF[k] /= strata_N[k];
                k++;
            }
		}

		// Print variant info
        for (int ii = 0; ii <= n_varInfo; ii++) {
        	oss << values[ii] << "\t";
        }

		if (strata) {
			for (int ii = 0; ii < n_cart; ii++) {
				oss << strata_N[ii] << "\t" << strata_AF[ii] << "\t";
			}
			std::fill(strata_N.begin(), strata_N.end(), 0);
			std::fill(strata_AF.begin(), strata_AF.end(), 0);			
		}

		oss << values[beta_marginal] << "\t";
		if (robust == 1) {
			oss << values[R_se_marginal] << "\t";
			if (printMeta || printFull) {
				oss << values[M_se_marginal] << "\t";
			}
		} else {
			oss << values[M_se_marginal] << "\t";
		}

		// Print betas
        for (int ii = printStart; ii < printEnd; ii++) {
        	oss << Beta_dot[ii] << "\t";
        }
		
		// Print robust variance-covariance first if applicable.
		if (robust == 1 ) {
			for (int ii = printStart; ii < printEnd; ii++) {
            	oss << sqrt(RV_dot[ii * Sq1 + ii]) << "\t";
        	}
        	for (int ii = printStart; ii < printEnd; ii++) {
            	for (int jj = printStart; jj < printEnd; jj++) {
                	if (ii < jj) {
                    	oss << RV_dot[ii * Sq1 + jj] << "\t";
                	}
           	 	}
        	}

			if (printMeta || printFull) {
				for (int ii = printStart; ii < printEnd; ii++) {
            		oss << sqrt(Vi_dot[ii * Sq1 + ii]) << "\t";
				}
				for (int ii = printStart; ii < printEnd; ii++) {
					for (int jj = printStart; jj < printEnd; jj++) {
						if (ii < jj) {
							oss << Vi_dot[ii * Sq1 + jj] << "\t";
						}
					}
				}
			}

			oss << values[R_pval_marginal] << "\t" <<  R_pvalInt << "\t" << R_pvalJoint;
			
			if (printMeta || printFull) {
				oss << "\t" << values[M_pval_marginal] << "\t" << M_pvalInt << "\t" << M_pvalJoint;
			}
			oss << "\n";

		} else {
			for (int ii = printStart; ii < printEnd; ii++) {
            	oss << sqrt(Vi_dot[ii * Sq1 + ii]) << "\t";
			}
			if (printMeta || printFull) {
				for (int ii = printStart; ii < printEnd; ii++) {
					for (int jj = printStart; jj < printEnd; jj++) {
						if (ii < jj) {
							oss << Vi_dot[ii * Sq1 + jj] << "\t";
						}
					}
				}
			}
			oss << values[M_pval_marginal] << "\t" << M_pvalInt << "\t" << M_pvalJoint << "\n";
		}
	}

	results << oss.str();
	oss.str(std::string());
	oss.clear();
	results.close();
	
	auto end_time = std::chrono::high_resolution_clock::now();
	printExecutionTime(start_time, end_time);
	return 0;
}


void printExecutionTime(std::chrono::high_resolution_clock::time_point start_time, std::chrono::high_resolution_clock::time_point end_time)
{

    auto execution_time_ms   = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    auto execution_time_sec  = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    auto execution_time_min  = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time).count();
    auto execution_time_hour = std::chrono::duration_cast<std::chrono::hours>(end_time - start_time).count();

    execution_time_ms   = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    execution_time_sec  = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    execution_time_min  = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time).count();
    execution_time_hour = std::chrono::duration_cast<std::chrono::hours>(end_time - start_time).count();

    cout << "Execution time... ";
    if (execution_time_hour > 0)
        cout << "" << execution_time_hour << "h, ";
    if (execution_time_min > 0)
        cout << "" << execution_time_min % 60 << "m, ";
    if (execution_time_sec > 0)
        cout << "" << execution_time_sec % 60 << "s, ";
    if (execution_time_ms > 0)
        cout << "" << execution_time_ms % long(1E+3) << " ms";
    cout << "\n";
    cout << "Done. \n";
    cout << "*********************************************************\n";
}

void printOutputHeader(int n_varInfo, int printStart, int printEnd, int robust, bool strata, vector<string> catE_headers, vector<string> covNames, string outStyle, string output) 
{
	std::ofstream results(output, std::ofstream::binary);

	bool printMeta = false;
	bool printFull = false;
    if (outStyle.compare("meta") == 0) {
        printMeta = true;
    } else if (outStyle.compare("full") == 0) {
        printFull = true;
    }
    results << "SNPID" << ((n_varInfo == 8) ? "\tRSID\t" : "\t") << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t";
    if (strata) {
        for (size_t i = 0; i < catE_headers.size(); i++) {
            results << "N_" << catE_headers[i] << "\t";
            results << "AF_" << catE_headers[i] << "\t";
        }
    }

    results << "Beta_Marginal" << "\t";
    if (robust == 1) {
        results << "robust_SE_Beta_Marginal" << "\t";
		if (printMeta || printFull) {
			 results << "SE_Beta_Marginal" << "\t";
		}
    } else {
		results << "SE_Beta_Marginal" << "\t";
	}

	// Print beta header
	for (int i = printStart; i < printEnd; i++) {
        results << "Beta_" << covNames[i] << "\t";
    }

	if (robust == 1) {
		for (int i = printStart; i < printEnd; i++) {
			results << "robust_SE_" << covNames[i] << "\t";  
		}

		for (int i = printStart; i < printEnd; i++) {
            for (int j = printStart; j < printEnd; j++) {
                if (i < j) {
                   results << "robust_Cov_Beta_" << covNames[i] << "_" << covNames[j] << "\t";  
                } 
            }
        }

		if (printMeta || printFull) {
			for (int i = printStart; i < printEnd; i++) {
                for (int j = printStart; j < printEnd; j++) {
                    if (i == j) {
                        results << "SE_Beta_" << covNames[j] << "\t"; 
                    }
                }
            }
            for (int i = printStart; i < printEnd; i++) {
                for (int j = printStart; j < printEnd; j++) {
                    if (i < j) {
                        results << "Cov_Beta_" << covNames[i] << "_" << covNames[j] << "\t"; 
                    }
                }
            }
		}

		results  << "robust_P_value_Marginal" << "\t" << "robust_P_Value_Interaction" << "\t" << "robust_P_Value_Joint" << "\t";
        (printMeta || printFull) ? results << "P_value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n" : results << "\n";
	} else {

		for (int i = printStart; i < printEnd; i++) {
            for (int j = printStart; j < printEnd; j++) {
                if (i == j) {
                    results << "SE_Beta_" << covNames[j] << "\t"; 
                }
            }
        }

		if (printMeta || printFull) {
            for (int i = printStart; i < printEnd; i++) {
                for (int j = printStart; j < printEnd; j++) {
                    if (i < j) {
                        results << "Cov_Beta_" << covNames[i] << "_" << covNames[j] << "\t"; 
                    }
                }
            }
		}
		results << "P_value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n";
	}

    results.close();
}

std::vector<std::string> cartesian_vec( vector<vector<string> >& v ) 
{
  auto product = []( long long a, vector<string>& b ) { return a*b.size(); };
  const long long N = accumulate( v.begin(), v.end(), 1LL, product );
  vector<string> u(v.size());
  std::vector<std::string> vec;
  for( long long n=0 ; n<N ; ++n ) {
    lldiv_t q { n, 0 };
    for( long long i=v.size()-1 ; 0<=i ; --i ) {
      q = div( q.quot, v[i].size() );
      u[i] = v[i][q.rem];
    }

    string tmp = "";
    for( size_t i = 0; i < u.size(); i++) { 
        tmp+=u[i]; 
    }
    vec.push_back(tmp);
  }
  return vec;
}

std::vector<std::string> cartesian_vec_sep( vector<vector<string> >& v) 
{
  auto product = []( long long a, vector<string>& b ) { return a*b.size(); };
  const long long N = accumulate( v.begin(), v.end(), 1LL, product );
  vector<string> u(v.size());
  std::vector<std::string> vec(N);
  for( long long n=0 ; n<N ; ++n ) {
    lldiv_t q { n, 0 };
    for( long long i=v.size()-1 ; 0<=i ; --i ) {
      q = div( q.quot, v[i].size() );
      u[i] = v[i][q.rem];
    }

    string tmp = "";
    for( size_t i = 0; i < u.size(); i++) { 
        tmp = (i == 0) ? u[i] : tmp + "_" + u[i]; 
    }
    vec[n] = tmp;
  }
  return vec;
}