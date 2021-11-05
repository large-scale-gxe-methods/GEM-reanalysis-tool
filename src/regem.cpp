#include "declars.h"

std::vector<std::string> cartesian_vec( vector<vector<string> >& v );
void printExecutionTime(std::chrono::high_resolution_clock::time_point start_time, std::chrono::high_resolution_clock::time_point end_time);
void printOutputHeader(int n_varInfo, int printStart, int printEnd, int robust, bool strata, bool printMeta, bool printFull, vector<string> catE_headers, vector<string> covNames, string output);
void compute_subcategorical(int Sq, std::vector<string> interactions, std::vector<string> stratum_names, bool* strata_ret, std::vector<int>* len_ret, std::vector<int>* idx_ret, std::vector<string>* header_ret);

int main(int argc, char* argv[]) 
{
	// Process command line arguments
	CommandLine cmd;
	cmd.processCommandLine(argc, argv);

	auto start_time = std::chrono::high_resolution_clock::now();

	int  printStart = cmd.printStart;
	int  printEnd   = cmd.printEnd;
	bool printMeta  = cmd.printMeta;
	bool printFull  = cmd.printFull;
	int  robust     = 0;

	int Sq     = cmd.numIntSelCol + cmd.numExpSelCol;
	int Sq1    = cmd.numIntSelCol + cmd.numExpSelCol + 1;
	int expSq  = cmd.numExpSelCol;
	int expSq1 = cmd.numExpSelCol + 1;

	std::vector<std::string> interactions = cmd.interactions;



	// Read input file
	std::ifstream res;
    	res.open(cmd.inFile);
    	if (!res.is_open()) {
        	cerr << "\nERROR: Cannot open results file. \n\n" << endl;
        	exit(1);
    	}

	// Get dispersion
	std::string line;
    	getline(res, line);

	double sigma2 = 0.0;
	if (line.rfind("#dispersion", 0) != 0) {
		cerr << "\nERROR: The first line in the input file is not '#dispersion'.\n\n";
		exit(1); 
	} else {
		std::istringstream iss(line);
		std::string value;
		iss >> value >> value;
		sigma2 = std::stod(value);
	}



	// Read the column headers
	int n_strata = 0;
	std::vector<std::string> stratum_names;
	std::unordered_map<std::string, int> column_names;

	getline(res, line);

	int header_i = 0;
	std::string header;
    	std::istringstream iss(line);
    	while (getline(iss, header, '\t')) {
        	header.erase(std::remove(header.begin(), header.end(), '\r'), header.end());
		if (column_names.find(header) != column_names.end()) {
            		cerr << "\nERROR: There are duplicate header names (" << header << ") in the results file.\n\n";
            		exit(1);
        	}
		column_names[header] = header_i;

        	if (header.rfind("Beta_G-", 0) == 0) {
            		string tmp_name = header;
            		tmp_name.erase(0, 7);
            		if (std::find(interactions.begin(), interactions.end(), tmp_name) == interactions.end()) {
                		interactions.push_back(tmp_name);
            		}
       		}

		if (header.rfind("N_", 0) == 0) {
			if (header != "N_Samples") {
				stratum_names.push_back(header);
				n_strata++;
			}
		}

		if ((robust !=1) && (header.rfind("robust_", 0) == 0)) {
			robust = 1;
		}

        	header_i++;
    	}


	// Compute new strata if needed
	size_t n_strata_col = n_strata * 2;
	bool strata = false;
	std::vector<int> subCategories_len;
	std::vector<int> subCategories_idx;
	std::vector<std::string> subCategories_header;
	if (stratum_names.size() > 0 ) {
		compute_subcategorical(Sq, interactions, stratum_names, &strata, &subCategories_len, &subCategories_idx, &subCategories_header);
	}
	size_t n_subStrata = subCategories_header.size();


	// Get column header indices
	for (size_t i = 0; i < interactions.size(); i++) {
		interactions[i] = "G-" + interactions[i];
	}
	interactions.insert(interactions.begin(), "G");

	int	n_varInfo    = column_names["AF"] + 1;
	int beta_marginal    = column_names["Beta_Marginal"];
	int mb_se_marginal   = column_names["SE_Beta_Marginal"];
	int rb_se_marginal   = column_names["robust_SE_Beta_Marginal"];
	int mb_pval_marginal = column_names["P_Value_Marginal"];
	int rb_pval_marginal = column_names["robust_P_Value_Marginal"];

	int dim = interactions.size();
	vector<int> beta_idx(dim);
	vector<int> mb_cov_idx(dim*dim);	
	vector<int> rb_cov_idx(dim*dim);

	for (int i = 0; i < dim; i++) {
		// Beta columns
		if (column_names.find("Beta_" + interactions[i]) == column_names.end()) {
			cerr << "\nERROR: Column header Beta_" + interactions[i] + " does not exist in input file.\n\n";
			exit(1);
		} else {
			beta_idx[i] = column_names["Beta_" + interactions[i]];
		}

		// Model-based SE columns
		if (column_names.find("SE_Beta_" + interactions[i]) == column_names.end()) {
			cerr << "\nERROR: Column header SE_Beta_" + interactions[i] + " does not exist in the input file.\n\n";
			exit(1);
		} else {
			mb_cov_idx[(i * dim) + i] = column_names["SE_Beta_" + interactions[i]];	
		}

		// Robust SE columns
		if (robust == 1) {
			if (column_names.find("robust_SE_Beta_" + interactions[i]) == column_names.end()) {
				cerr << "\nERROR: Column header robust_SE_Beta_" + interactions[i] + " does not exist in the input file.\n\n";
				exit(1);
			} else {
				rb_cov_idx[(i * dim) + i] = column_names["robust_SE_Beta_" + interactions[i]];	
			}
		}

		// Model-based Covariance columns
		for (int j = (i+1); j < dim; j++) {
			if (column_names.find("Cov_Beta_" + interactions[i] + "_" + interactions[j]) != column_names.end()) {
				int tmp_idx = column_names["Cov_Beta_" + interactions[i] + "_" + interactions[j]];
				mb_cov_idx[(i * dim) + j] = tmp_idx;
				mb_cov_idx[(j * dim) + i] = tmp_idx;
			} 
			else if (column_names.find("Cov_Beta_" + interactions[j] + "_" + interactions[i]) != column_names.end()) {
 				int tmp_idx = column_names["Cov_Beta_" + interactions[j] + "_" + interactions[i]];
				mb_cov_idx[(i * dim) + j] = tmp_idx;
				mb_cov_idx[(j * dim) + i] = tmp_idx;
			} else {
				cerr << "\nERROR: Column header Cov_Beta_" + interactions[i] + "_" + interactions[j] + " does not exist in the input file.\n\n";
				exit(1);
			}
		}

		// Robust Covariance columns
		if (robust == 1) {
			for (int j = (i+1); j < dim; j++) {
				if (column_names.find("robust_Cov_Beta_" + interactions[i] + "_" + interactions[j]) != column_names.end()) {
					int tmp_idx = column_names["robust_Cov_Beta_" + interactions[i] + "_" + interactions[j]];
					rb_cov_idx[(i * dim) + j] = tmp_idx;
					rb_cov_idx[(j * dim) + i] = tmp_idx;
				}
				else if (column_names.find("robust_Cov_Beta_" + interactions[j] + "_" + interactions[i]) != column_names.end()) {
					int tmp_idx = column_names["robust_Cov_Beta_" + interactions[j] + "_" + interactions[i]];
					rb_cov_idx[(i * dim) + j] = tmp_idx;
					rb_cov_idx[(j * dim) + i] = tmp_idx;					
				}
				else {
					cerr << "\nERROR: Column header Cov_Beta_" + interactions[i] + "_" + interactions[j] + " does not exist in the input file.\n\n";
					exit(1);
				}
			}	
		}
	}



	// Check if there are results in the input file
    	int nvars = 0;
    	while (getline(res, line)) nvars++;
	if (nvars == 0) {
		cerr << "\nERROR: No GEM results in the input file.\n\n.";
		exit(1);
	} else {
		res.clear(); 
		res.seekg(0, res.beg);
    		getline(res, line);
		getline(res, line);	
	}


	// Create the output file and write column header names
	std::ofstream results(cmd.outFile, std::ios_base::app);
	std::ostringstream oss;
	printOutputHeader(n_varInfo, cmd.printStart, cmd.printEnd, robust, strata, cmd.printMeta, cmd.printFull, subCategories_header, interactions, cmd.outFile);





	// REGEM begins
	vector<double> beta_vec(dim);
	vector<double> mb_vec(dim*dim);
	vector<double> rb_vec(dim*dim);
	vector<double> strata_N(n_subStrata);
	vector<double> strata_AF(n_subStrata);
	double* Beta = &beta_vec[0];
	double* MV   = &mb_vec[0];
	double* RV   = &rb_vec[0];

	double* StempE   = new double[expSq];
	double* StempGE  = new double[expSq1];
	double* U        = new double[dim];
	double* Beta_dot = new double[Sq1];
	double* V        = new double[dim * dim];
	double* Vi_dot   = new double[Sq1*Sq1];
	double* VE_dot   = new double[expSq*expSq];
	double* O        = new double[dim*dim];
	double* O_dot    = new double[Sq1*Sq1];
	double* RV_dot   = new double[Sq1*Sq1];
	double* A_i      = new double[expSq1 * expSq1];
	double M_pvalM, M_pvalInt, M_pvalJoint, R_pvalM, R_pvalInt, R_pvalJoint;

	boost::math::chi_squared chisq_dist_M(1);
	boost::math::chi_squared chisq_dist_Int(expSq);
	boost::math::chi_squared chisq_dist_Joint(expSq1);

	int cnt = 0;
	for (int p = 0; p < nvars; p++) {

		getline(res, line);
		std::istringstream iss(line);
		string value;
		vector <string> values;
		while (getline(iss, value, '\t')) values.push_back(value);

		for (int i = 0; i < dim; i++) {
			sscanf(values[beta_idx[i]].c_str(), "%lf", &beta_vec[i]);
			for (int j = 0; j < dim; j++) {
	 			sscanf(values[mb_cov_idx[(i * dim) + j]].c_str(), "%lf", &mb_vec[(i * dim) + j]);
	 		}
			mb_vec[(i * dim) + i] *= mb_vec[(i * dim) + i];
		}

		if (robust == 1) {
			for (int i = 0; i < dim; i++) {
				for (int j = 0; j < dim; j++) {
					sscanf(values[rb_cov_idx[(i * dim) + j]].c_str(), "%lf", &rb_vec[(i * dim) + j]);
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

		if (strata) 
		{
			size_t i = 0;
			size_t k = 0;
			
			vector<double> tv(n_strata_col);
			for (int j = 0; j < (n_strata_col); j++) 
				tv[j] = std::stod(values[(n_varInfo) + j]);

            		while (i < subCategories_idx.size()) 
			{
                		for (int j = 0; j < subCategories_len[k]; j++) 
				{
                    			strata_N[k]  += tv[subCategories_idx[i]];
                    			strata_AF[k] += (tv[subCategories_idx[i]] * ((!std::isnan(tv[subCategories_idx[i]+1])) ? tv[subCategories_idx[i]+1] : 0));
                    			i++;
                		}
                		strata_AF[k] = (strata_N[k] != 0) ? strata_AF[k] /= strata_N[k] : NAN;
                		k++;
            		}
		}



		// Print variant info
        	for (int ii = 0; ii < n_varInfo; ii++) {
        		oss << values[ii] << "\t";
        	}

		if (strata) {
			for (int ii = 0; ii < n_subStrata; ii++) {
				oss << strata_N[ii] << "\t" << strata_AF[ii] << "\t";
			}
			std::fill(strata_N.begin(), strata_N.end(), 0);
			std::fill(strata_AF.begin(), strata_AF.end(), 0);			
		}

		oss << values[beta_marginal] << "\t";
		if (robust == 1) {
			oss << values[rb_se_marginal ] << "\t";
			if (printMeta || printFull) {
				oss << values[mb_se_marginal] << "\t";
			}
		} else {
			oss << values[mb_se_marginal] << "\t";
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

			oss << values[rb_pval_marginal ] << "\t" <<  R_pvalInt << "\t" << R_pvalJoint;
			
			if (printMeta || printFull) {
				oss << "\t" << values[mb_pval_marginal ] << "\t" << M_pvalInt << "\t" << M_pvalJoint;
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
			oss << values[mb_pval_marginal] << "\t" << M_pvalInt << "\t" << M_pvalJoint << "\n";
		}
		cnt++;

		if (cnt % 100000 != 0) {
        	results << oss.str();
        	oss.str(std::string());
        	oss.clear();
			
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



void printOutputHeader(int n_varInfo, int printStart, int printEnd, int robust, bool strata, bool printMeta, bool printFull, vector<string> catE_headers, vector<string> covNames, string output) 
{
	std::ofstream results(output, std::ofstream::binary);

    	results << "SNPID" << ((n_varInfo == 8) ? "\tRSID\t" : "\t") << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t";
    	if (strata) {
        	for (size_t i = 0; i < catE_headers.size(); i++) {
            		results << "N_"  << catE_headers[i] << "\t";
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

		if (printMeta || printFull) 
		{
            		for (int i = printStart; i < printEnd; i++) 
			{
                		for (int j = printStart; j < printEnd; j++) 
				{
                    			if (i < j) 
					{
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
	
	for (long long n=0 ; n<N ; ++n ) 
  	{
		lldiv_t q { n, 0 };
    	for (long long i=v.size()-1 ; 0<=i ; --i )
		{
    		q = lldiv( q.quot, v[i].size() );
    		u[i] = v[i][q.rem];
    	}

    	string tmp = "";
    	for (size_t i = 0; i < u.size(); i++) 
		{ 
        	tmp = (i == 0) ? u[i] : tmp + "_" + u[i]; 
    	}
    	vec.push_back(tmp);
  	}
  	return vec;
}


void compute_subcategorical(int Sq, std::vector<string> interactions, std::vector<string> stratum_names, bool* strata_ret, std::vector<int>* len_ret, std::vector<int>* idx_ret, std::vector<string>* header_ret)
{
	std::vector<std::string> categorical_names;
	std::vector<std::string> header_split;

	// Get the first stratum header name in the results file.
	std::string header;
	std::stringstream ss(stratum_names[0]);
	while(std::getline(ss, header, '_')) 
	{
		header_split.push_back(header);
	}
	size_t split_len = header_split.size();

	// Look for the interactions in the name of the stratum header name (Slow process).
    	size_t k = 1;
	size_t j = split_len;
	int n_categories = 0;
	while (k < j) 
	{
		int m = 1;
        	bool is_cat = false;
		while (!is_cat) 
		{
			int counter = 0;
            		std::string s = header_split[k];
            		for (int i = k+1; i < header_split.size() - m; i++) 
			{
                		s += "_" + header_split[i];
                		counter++;
            		}
            		if (std::find(interactions.begin(), interactions.end(), s) != interactions.end()) 
			{
				is_cat = true;
                		categorical_names.push_back(s);
			    	k = k + counter + 1;
			    	j--;
				n_categories++;
            		}	
            		m++;
        	}
	}
	
	// Get indices of the sub categorical interaction terms in results file
	int n_subCategories = 0;
	std::string subcat_string;
	std::vector<int> subCategories_idx;
	std::vector<std::string> subCategories;
	for (int i = 0; i < Sq; i++) 
	{
		auto it = std::find(categorical_names.begin(), categorical_names.end(), interactions[i]);
		if (it != categorical_names.end()) 
		{
			int index = it - categorical_names.begin();
			subCategories_idx.push_back(split_len - (n_categories - index));
			subCategories.push_back(interactions[i]);
			subcat_string += interactions[i] + "_";
			n_subCategories++;
		}
	}

	if (n_subCategories > 0 ) 
	{
		*strata_ret = true;

		// Get each sub categorical interaction term's unique categories based on the header name
		std::vector<std::vector<std::string>> tmp_vec(stratum_names.size());
		std::vector<std::vector<std::string>> subCategories_strata(n_subCategories);
		for (size_t i = 0; i < stratum_names.size(); i++) 
		{
			std::string line2;
			std::stringstream ss(stratum_names[i]);
			while(std::getline(ss, line2, '_')) 
			{
				tmp_vec[i].push_back(line2);
			}

			for (int j = 0; j < n_subCategories; j++) 
			{
				if (std::find(subCategories_strata[j].begin(), subCategories_strata[j].end(), tmp_vec[i][subCategories_idx[j]]) == subCategories_strata[j].end()) 
				{
					subCategories_strata[j].push_back(tmp_vec[i][subCategories_idx[j]]);
				}
			}

		}

		// Compute the cartesian product between the sub categorical interaction terms
		std::vector<std::string> subCategories_cartesian = cartesian_vec(subCategories_strata);
	
		std::vector<int> subCategories_Nidx;
		std::vector<int> subCategories_count;
		for(size_t i = 0; i < subCategories_cartesian.size(); i++) 
		{
			int count = 0;
			for (size_t j = 0; j < stratum_names.size(); j++)
			{
				string tmp = "";
				for (int k = 0; k < n_subCategories; k++) 
				{
					tmp = (k == 0) ? tmp_vec[j][subCategories_idx[k]] : tmp + "_" + tmp_vec[j][subCategories_idx[k]];
				}
				if (tmp == subCategories_cartesian[i]) 
				{
					subCategories_Nidx.push_back(j*2);
					count++;
				} 
			}
			subCategories_count.push_back(count);
		}

		vector<string> subCategories_header;
		for (size_t i = 0; i < subCategories_cartesian.size(); i++) {
			subCategories_header.push_back(subcat_string + subCategories_cartesian[i]);
		}

		*header_ret = subCategories_header;
		*idx_ret = subCategories_Nidx;
		*len_ret = subCategories_count;
	}	

}