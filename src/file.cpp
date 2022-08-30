#include "regem.h"

void processFileHeader(std::string fileName, std::vector<std::string> interactions, size_t nInt1, FileInfo* fip) 
{
    // Read input file
    std::ifstream file;
    file.open(fileName);
    if (!file.is_open()) {
        cout << "\nERROR: Cannot open input file.\n\n" << endl;
        exit(1);
    }

    // Ignore comments in the text file except dispersion line
    std::string line;
    std::string header;
    while(getline(file, line))
    {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.rfind("#", 0) != 0) {
            break; 
        }

        if (line.rfind("#dispersion:", 0) == 0) {
            std::istringstream iss(line);
            std::string value;
            iss >> value >> value;
            fip->sigma2 = std::stod(value);
        }
    }


    // Read the column headers
    std::vector<std::string> categorical_columns;
    std::unordered_map<std::string, int> columnNames;

    int header_i = 0;
    std::istringstream iss(line);
    while (getline(iss, header, '\t')) {
        if (columnNames.find(header) != columnNames.end()) {
            cout << "\nERROR: There are duplicate header names (" << header << ") in the results file.\n\n";
            exit(1);
        }
        columnNames[header] = header_i;
        
        if (header.rfind("Beta_G-", 0) == 0) {
            header.erase(0, 7);
            if (std::find(interactions.begin(), interactions.end(), header) == interactions.end()) {
                interactions.push_back(header);
            }
        }

        if ((header.rfind("N_", 0) == 0) && (header != "N_Samples")) {
            header.erase(0, 2);
            categorical_columns.push_back(header);
        }

        if ((!fip->robust) && (header.rfind("robust_Cov", 0) == 0)) {
            fip->robust = true;
        }
        header_i++;
    }


    // Get variant info indices
    get_variantInfo_indices(columnNames, fip);

    // Generate new headers for the categorical variables
    generate_subcategorical_columns(nInt1, columnNames, interactions, categorical_columns, fip);

    // Get marginal summary statistics indices
    get_marginal_indices(columnNames, fip);


    // Get indices of interaction terms
    size_t dim = interactions.size() + 1;
    fip->betaIntIndex.resize(dim);
    fip->mbCovIndex.resize(dim * dim);
    fip->rbCovIndex.resize(dim * dim);

    interactions.insert(interactions.begin(), "G");
    for (size_t i = 1; i < dim; i++) {
        interactions[i] = "G-" + interactions[i];
    }

    for (size_t i = 0; i < dim; i++) {
        // Beta columns
        if (columnNames.find("Beta_" + interactions[i]) == columnNames.end()) {
            cout << "\nERROR: Column header Beta_" + interactions[i] + " does not exist in input file.\n\n";
            exit(1);
        } else {
            fip->betaIntIndex[i] = columnNames["Beta_" + interactions[i]];
        }

        // Model-based SE columns
        if (columnNames.find("SE_Beta_" + interactions[i]) == columnNames.end()) {
            cout << "\nERROR: Column header SE_Beta_" + interactions[i] + " does not exist in the input file.\n\n";
            exit(1);
        } else {
            fip->mbCovIndex[(i * dim) + i] = columnNames["SE_Beta_" + interactions[i]];	
        }

        // Robust SE columns
        if (fip->robust) {
            if (columnNames.find("robust_SE_Beta_" + interactions[i]) == columnNames.end()) {
                cout << "\nERROR: Column header robust_SE_Beta_" + interactions[i] + " does not exist in the input file.\n\n";
                exit(1);
            } else {
                fip->rbCovIndex[(i * dim) + i] = columnNames["robust_SE_Beta_" + interactions[i]];	
            }
        }

        // Model-based Covariance columns
        for (size_t j = (i+1); j < dim; j++) {
            if (columnNames.find("Cov_Beta_" + interactions[i] + "_" + interactions[j]) != columnNames.end()) {
                int tmp_idx = columnNames["Cov_Beta_" + interactions[i] + "_" + interactions[j]];
                fip->mbCovIndex[(i * dim) + j] = tmp_idx;
                fip->mbCovIndex[(j * dim) + i] = tmp_idx;
            } 
            else if (columnNames.find("Cov_Beta_" + interactions[j] + "_" + interactions[i]) != columnNames.end()) {
                int tmp_idx = columnNames["Cov_Beta_" + interactions[j] + "_" + interactions[i]];
                fip->mbCovIndex[(i * dim) + j] = tmp_idx;
                fip->mbCovIndex[(j * dim) + i] = tmp_idx;
            } else {
                cout << "\nERROR: Column header Cov_Beta_" + interactions[i] + "_" + interactions[j] + " does not exist in the input file.\n\n";
                exit(1);
            }
        }

        // Robust Covariance columns
        if (fip->robust) {
            for (size_t j = (i+1); j < dim; j++) {
                if (columnNames.find("robust_Cov_Beta_" + interactions[i] + "_" + interactions[j]) != columnNames.end()) {
                    int tmp_idx = columnNames["robust_Cov_Beta_" + interactions[i] + "_" + interactions[j]];
                    fip->rbCovIndex[(i * dim) + j] = tmp_idx;
                    fip->rbCovIndex[(j * dim) + i] = tmp_idx;
                }
                else if (columnNames.find("robust_Cov_Beta_" + interactions[j] + "_" + interactions[i]) != columnNames.end()) {
                    int tmp_idx = columnNames["robust_Cov_Beta_" + interactions[j] + "_" + interactions[i]];
                    fip->rbCovIndex[(i * dim) + j] = tmp_idx;
                    fip->rbCovIndex[(j * dim) + i] = tmp_idx;					
                }
                else {
                    cout << "\nERROR: Column header Cov_Beta_" + interactions[i] + "_" + interactions[j] + " does not exist in the input file.\n\n";
                    exit(1);
                }
            }
        }
    }
    fip->dim = dim;
    fip->interaction_names = interactions;

    file.close();
}


void get_variantInfo_indices(std::unordered_map<std::string, int> columnNames, FileInfo* fip)
{
    if (columnNames.find("SNPID") != columnNames.end()) 
    {
        fip->varInfoIndices.push_back(columnNames["SNPID"]);
    } 
    else 
    {
        cout << "\nERROR: SNPID column does not exists.\n\n";
        exit(1);
    }

    if (columnNames.find("RSID") != columnNames.end()) 
    {
        fip->varInfoIndices.push_back(columnNames["RSID"]);
    }

    if (columnNames.find("CHR") != columnNames.end())
    {
        fip->varInfoIndices.push_back(columnNames["CHR"]);
    }
    else
    {
        cout << "\nERROR: CHR column does not exists.\n\n";
        exit(1);
    }

    if (columnNames.find("POS") != columnNames.end())
    {
        fip->varInfoIndices.push_back(columnNames["POS"]);
    }
    else
    {
        cout << "\nERROR: POS column does not exists.\n\n";
        exit(1);
    }

    if (columnNames.find("Non_Effect_Allele") != columnNames.end())
    {
        fip->varInfoIndices.push_back(columnNames["Non_Effect_Allele"]);
    }
    else
    {
        cout << "\nERROR: Non_Effect_Allele column does not exists.\n\n";
        exit(1);
    }

    if (columnNames.find("Effect_Allele") != columnNames.end())
    {
        fip->varInfoIndices.push_back(columnNames["Effect_Allele"]);
    }
    else
    {
        cout << "\nERROR: Effect_Allele column does not exists.\n\n";
        exit(1);
    }

    if (columnNames.find("N_Samples") != columnNames.end())
    {
        fip->varInfoIndices.push_back(columnNames["N_Samples"]);
    }
    else
    {
        cout << "\nERROR: N_Samples column does not exists.\n\n";
        exit(1);
    }

    if (columnNames.find("AF") != columnNames.end())
    {
        fip->varInfoIndices.push_back(columnNames["AF"]);
    }
    else
    {
        cout << "\nERROR: AF column does not exists.\n\n";
        exit(1);
    }
    fip->nVarInfo = fip->varInfoIndices.size();
}

void get_marginal_indices(std::unordered_map<std::string, int> columnNames, FileInfo* fip) 
{
    // Get marginal summary statistics indices
    if (columnNames.find("Beta_Marginal") != columnNames.end())
    {
        fip->betaMarginalColumn = columnNames["Beta_Marginal"];
    }
    else
    {
        cout << "\nERROR: Beta_Marginal column does not exists.\n\n";
        exit(1);
    } 

    if (columnNames.find("SE_Beta_Marginal") != columnNames.end())
    {
        fip->mbSeMarginalColumn = columnNames["SE_Beta_Marginal"];
    }
    else
    {
        cout << "\nERROR: SE_Beta_Marginal column does not exists.\n\n";
        exit(1);
    } 

    if (columnNames.find("P_Value_Marginal") != columnNames.end())
    {
        fip->mbPvalMarginalColumn = columnNames["P_Value_Marginal"];
    }
    else
    {
        cout << "\nERROR: P_Value_Marginal column does not exists.\n\n";
        exit(1);
    } 

    if (fip->robust) {
        if (columnNames.find("robust_SE_Beta_Marginal") != columnNames.end())
        {
            fip->rbSeMarginalColumn = columnNames["robust_SE_Beta_Marginal"];
        }
        else
        {
            cout << "\nERROR: robust_SE_Beta_Marginal column does not exists.\n\n";
            exit(1);
        } 

        if (columnNames.find("robust_P_Value_Marginal") != columnNames.end())
        {
            fip->rbPvalMarginalColumn = columnNames["robust_P_Value_Marginal"];
        }
        else
        {
            cout << "\nERROR: robust_P_Value_Marginal column does not exists.\n\n";
            exit(1);
        } 
    }
}


std::vector<std::string> cartesian_product(std::vector<std::vector<std::string>>& v ) 
{
    auto product = []( long long a, std::vector<std::string>& b ) { return a*b.size(); };
    const long long N = std::accumulate( v.begin(), v.end(), 1LL, product );

    std::vector<std::string> u(v.size());
    std::vector<std::string> vec;
    for (long long n=0 ; n<N ; ++n ) {
        lldiv_t q { n, 0 };
        for (long long i=v.size()-1 ; 0<=i ; --i ) {
            q = lldiv( q.quot, v[i].size() );
            u[i] = v[i][q.rem];
        }

        std::string temp = "";
        for (size_t i = 0; i < u.size(); i++) { 
            temp = (i == 0) ? u[i] : temp + "_" + u[i]; 
        }
        vec.push_back(temp);
    }
    return vec;
}


void generate_subcategorical_columns(size_t nInt1, std::unordered_map<std::string, int> columnNames, std::vector<std::string> interactions, std::vector<std::string> categorical_columns, FileInfo* fip)
{
    // Get the first stratum header name in the results file.
    std::string header;
    std::vector<std::string> header_split;
    std::stringstream ss(categorical_columns[0]);
    while(std::getline(ss, header, '_')) 
    {
        header_split.push_back(header);
    }
    size_t split_length = header_split.size();


    // Get the count of all categorical variables
    std::vector<std::string> categorical_names;
    for (size_t i = 0; i < split_length; i++) {
        if (std::find(interactions.begin(), interactions.end(), header_split[i]) != interactions.end()) {
            categorical_names.push_back(header_split[i]);
        }
    }
    size_t n = categorical_names.size();

    // Check if the selected exposure or int covar are categorical
    std::string subcategorical_prefix = "";
    std::vector<size_t> subcategorical_index;
    std::vector<std::string> subcategorical_names;
    for (size_t i = 0; i < (nInt1 - 1); i++) {
        auto it = std::find(categorical_names.begin(), categorical_names.end(), interactions[i]);
        if (it != categorical_names.end()) {
            int index = it - categorical_names.begin();
            subcategorical_index.push_back(split_length - (n - index));
            subcategorical_names.push_back(interactions[i]);
            subcategorical_prefix += interactions[i] + "_";
        }
    }
    size_t nsubcategorical = subcategorical_names.size();


    if (nsubcategorical > 0 ) {
        fip->subcategorical_exists = true;

        // Get unique strata of all subcategorical variables based on the header name
        std::vector<std::vector<std::string>> subcategorical_strata(nsubcategorical);

        std::vector<std::vector<std::string>> temp(categorical_columns.size());
        for (size_t i = 0; i < categorical_columns.size(); i++) {
            std::stringstream ss(categorical_columns[i]);
            while(std::getline(ss, header, '_')) 
            {
                temp[i].push_back(header);
            }

            for (size_t j = 0; j < nsubcategorical; j++) {
                if (std::find(subcategorical_strata[j].begin(), subcategorical_strata[j].end(), temp[i][subcategorical_index[j]]) == subcategorical_strata[j].end()) {
                    subcategorical_strata[j].push_back(temp[i][subcategorical_index[j]]);
                }
            }
        }

        // Compute the cartesian product between the sub categorical interaction terms
        std::vector<std::string> subcategorical_cp = cartesian_product(subcategorical_strata);

        fip->subcategorical_file_index.resize(subcategorical_cp.size());
        for(size_t i = 0; i < subcategorical_cp.size(); i++) {
            for (size_t j = 0; j < categorical_columns.size(); j++) {
                std::string cp = temp[j][subcategorical_index[0]];
                for (size_t k = 1; k < nsubcategorical; k++) {
                    cp += "_" + temp[j][subcategorical_index[k]];
                }

                if (cp == subcategorical_cp[i]) {
                    fip->subcategorical_file_index[i].push_back(columnNames["N_" + categorical_columns[j]]);
                    fip->subcategorical_file_index[i].push_back(columnNames["AF_" + categorical_columns[j]]);
                } 
            }
        }

        for (size_t i = 0; i < subcategorical_cp.size(); i++) {
            fip->subcategorical_header.push_back(subcategorical_prefix + subcategorical_cp[i]);
        }
    }
}


void printOutputHeader(FileInfo* fip,  CommandLine* cmd) 
{
    size_t printStart = cmd->printStart;
    size_t printEnd   = cmd->printEnd;
    std::vector<std::string> interaction_names = fip->interaction_names;
    std::ofstream results(cmd->outFile, std::ofstream::binary);

    if (cmd->outStyle.compare("full")) {
        results << "#dispersion: " << fip->sigma2 << "\n";
    }
    results << "SNPID" << ((fip->nVarInfo == 8) ? "\tRSID\t" : "\t") << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t";
    
    if (fip->subcategorical_exists) {
        std::vector<std::string> subcategorical_header = fip->subcategorical_header;
        for (size_t i = 0; i < subcategorical_header.size(); i++) {
            results << "N_"  << subcategorical_header[i] << "\t";
            results << "AF_" << subcategorical_header[i] << "\t";
        }
    }
    
    results << "Beta_Marginal" << "\t";
    if (fip->robust) {
        results << "SE_Beta_Marginal" << "\t";
        results << "robust_SE_Beta_Marginal" << "\t";
    } else {
        results << "SE_Beta_Marginal" << "\t";
    }


    // Print Int beta header
    for (size_t i = printStart; i < printEnd; i++) {
        results << "Beta_" << interaction_names[i] << "\t";
    }

    if (fip->robust) {
        for (size_t i = printStart; i < printEnd; i++) {
            results << "SE_Beta_" << interaction_names[i] << "\t";
        }
        for (size_t i = printStart; i < printEnd; i++) {
            for (size_t j = printStart; j < printEnd; j++) {
                if (i < j) {
                    results << "Cov_Beta_" << interaction_names[i] << "_" << interaction_names[j] << "\t"; 
                }
            }
        }

        for (size_t i = printStart; i < printEnd; i++) {
            results << "robust_SE_Beta_" << interaction_names[i] << "\t";  
        }
        for (size_t i = printStart; i < printEnd; i++) {
            for (size_t j = printStart; j < printEnd; j++) {
                if (i < j) {
                    results << "robust_Cov_Beta_" << interaction_names[i] << "_" << interaction_names[j] << "\t";  
                } 
            }
        }

        results << "P_Value_Marginal\t" << "P_Value_Interaction\t" << "P_Value_Joint\t";
        results  << "robust_P_Value_Marginal\t" << "robust_P_Value_Interaction\t" << "robust_P_Value_Joint\n";

    } else {
        for (size_t i = printStart; i < printEnd; i++) {
            results << "SE_Beta_" << interaction_names[i] << "\t";
        }

        for (size_t i = printStart; i < printEnd; i++) {
            for (size_t j = printStart; j < printEnd; j++) {
                if (i < j) {
                    results << "Cov_Beta_" << interaction_names[i] << "_" << interaction_names[j] << "\t"; 
                }
            }
        }
        results << "P_value_Marginal\t" << "P_Value_Interaction\t" << "P_Value_Joint\n";
    }
    results.close();
}