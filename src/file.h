#ifndef FILE_H
#define FILE_H

class FileInfo {
    public:
        std::vector<std::string> stratum_names;

        int dim;
        double sigma2;

        bool robust = false;
        size_t nVarInfo = 0;
        int betaMarginalColumn;
        int mbSeMarginalColumn;
        int rbSeMarginalColumn;
        int mbPvalMarginalColumn;
        int rbPvalMarginalColumn;
        std::vector<int> betaIntIndex;
        std::vector<int> mbCovIndex;
        std::vector<int> rbCovIndex;
        std::vector<int> varInfoIndices;
        std::vector<std::string> interaction_names;

        bool subcategorical_exists = false;
        std::vector<std::string> subcategorical_header;
        std::vector<std::vector<int>> subcategorical_file_index;
};

std::vector<std::string> cartesian_product(std::vector<std::vector<std::string>>& v );

void processFileHeader(std::string fileName, std::vector<std::string> interactions, size_t nInt1, FileInfo* fip);

void printOutputHeader(FileInfo* fip,  CommandLine* cmd); 

void get_variantInfo_indices(std::unordered_map<std::string, int> columnNames, FileInfo* fip);

void get_marginal_indices(std::unordered_map<std::string, int> columnNames, FileInfo* fip);

void generate_subcategorical_columns(size_t nInt1, std::unordered_map<std::string, int> columnNames, std::vector<std::string> interactions, std::vector<std::string> categorical_columns, FileInfo* fip);

#endif