// <REGEM: RE-analysis of GEM summary statistics>
// Copyright (C) <2021-2023> Duy T. Pham and Han Chen  

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.



#include "regem.h"
#include "time.h"

int main(int argc, char* argv[]) 
{
    // Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

    // Process command line arguments
    CommandLine* cmd = new CommandLine();
    cmd->processCommandLine(argc, argv);

    regem(cmd);

    // Stop timers
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    printTimeCompleted(wall0, wall1, cpu0, cpu1);
    
    return(0);
}

void regem(CommandLine* cmd) 
{
    // Process header file
    FileInfo* fip = new FileInfo();
    processFileHeader(cmd->inFile, cmd->interactions, cmd->nInt1, fip);


    // Read input file
    std::string line;
    std::string value;
    std::ifstream file;
    file.open(cmd->inFile);
    if (!file.is_open()) {
        cout << "\nERROR: Cannot open input file.\n\n" << endl;
        exit(1);
    }

    while(getline(file, line))
    {
        if (line.rfind("#", 0) != 0) {
            break;
        }
    }

    // Create the output file and write column header names
    printOutputHeader(fip, cmd);
    std::ofstream results(cmd->outFile, std::ios_base::app);
    std::ostringstream oss;




    // REGEM begins
    cout << "Running REGEM...\n";

    bool robust   = fip->robust;
    double sigma2 = fip->sigma2;
    size_t nExp   = cmd->nExp;
    size_t nExp1  = cmd->nExp1;
    size_t nInt1  = cmd->nInt1;
    size_t printStart = cmd->printStart;
    size_t printEnd   = cmd->printEnd;
    int centerIn   = cmd->centerIn;
    int centerOut   = cmd->centerOut;
    std::unordered_map<std::string, double> meanValues = cmd->meanValues;
	
    bool subcategorical_exists = fip->subcategorical_exists;
    size_t nsubcategorical_columns = fip->subcategorical_file_index.size();
    std::vector<std::vector<int>> subcategorical_file_index = fip->subcategorical_file_index;
    std::vector<double> strata_N(nsubcategorical_columns);
    std::vector<double> strata_AF(nsubcategorical_columns); 

    size_t dim = fip->dim;
    size_t nVarInfo = fip->nVarInfo;
    int betaMarginalColumn   = fip->betaMarginalColumn;
    int mbSeMarginalColumn   = fip->mbSeMarginalColumn;
    int rbSeMarginalColumn   = fip->rbSeMarginalColumn;
    int mbPvalMarginalColumn = fip->mbPvalMarginalColumn;
    int rbPvalMarginalColumn = fip->rbPvalMarginalColumn;
    std::vector<int> mbCovIndex     = fip->mbCovIndex;
    std::vector<int> rbCovIndex     = fip->rbCovIndex;
    std::vector<int> betaIntIndex   = fip->betaIntIndex;
    std::vector<int> varInfoIndices = fip->varInfoIndices;
    std::vector<std::string> interaction_names = fip->interaction_names;

    std::vector<double> beta(dim);
    std::vector<double> mb_v(dim*dim);
    std::vector<double> rb_v(dim*dim);
    double* Vi_dot   = new double[nInt1*nInt1];
    double* VE_dot   = new double[nExp*nExp];
    double* O        = new double[dim*dim];
    double* O_dot    = new double[nInt1*nInt1];
    double* MRV_prod = new double[dim*dim];
    double* OV_prod  = new double[nInt1 * nInt1];
    double* RV_dot   = new double[nInt1*nInt1];
    double* A_i      = new double[nExp1 * nExp1];
    double M_pvalInt, M_pvalJoint, R_pvalInt, R_pvalJoint;

    boost::math::chi_squared chisq_dist_M(1);
    boost::math::chi_squared chisq_dist_Int(nExp);
    boost::math::chi_squared chisq_dist_Joint(nExp1);

    int p = 0;
    while (getline(file, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        std::istringstream iss(line);
        std::vector<std::string> values;
        while (getline(iss, value, '\t')) values.push_back(value);

        // Betas
        for (size_t i = 0; i < dim; i++) {
            beta[i] = std::stod(values[betaIntIndex[i]]);
        }

        // MB covariances
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                mb_v[(i * dim) + j] = std::stod(values[mbCovIndex[(i * dim) + j]]);
            }
            mb_v[(i * dim) + i] *= mb_v[(i * dim) + i];
        }

        // Robust covariances
        if (robust) {
            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {
                    rb_v[(i * dim) + j] = std::stod(values[rbCovIndex[(i * dim) + j]]);
                }
                rb_v[(i * dim) + i] *= rb_v[(i * dim) + i];
            }
        }

	// Centering conversions of Betas and covariances
	if(centerIn != centerOut){
	    centerConversion(beta, mb_v, rb_v, robust, meanValues, interaction_names, dim, nExp, centerIn, centerOut);
	}
   
        // Compute V from model-based variance-covariance matrix
        matInv(&mb_v[0], dim);
        for (size_t i = 0; i < dim*dim; i++)
            mb_v[i] *= sigma2;

        // Compute U
        std::vector<double> u(dim, 0.0);
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                u[j] += (mb_v[(dim * j) + k] * beta[k]);
            }
        }

        subMatrix(&mb_v[0], Vi_dot, nInt1, nInt1, dim, nInt1, 0);
        matInv(Vi_dot, nInt1);

        // Compute betas
        std::vector<double> beta_dot(dim, 0.0);
        for (size_t j = 0; j < nInt1; j++) {
            for (size_t k = 0; k < nInt1; k++) {
                beta_dot[j] += (Vi_dot[(nInt1 * j) + k] * u[k]);
            }
        }


        if (robust) {
            //Compute Omega
            matNmatNprod(&mb_v[0], &rb_v[0], MRV_prod, dim, dim, dim);
            matNmatNprod(MRV_prod, &mb_v[0], O, dim, dim, dim);
            subMatrix(O, O_dot, nInt1, nInt1, dim, nInt1, 0);

            // Compute robust variance-covariance
            matNmatNprod(O_dot, Vi_dot,  OV_prod, nInt1, nInt1, nInt1);
            matNmatNprod(Vi_dot, OV_prod, RV_dot, nInt1, nInt1, nInt1);

            // Robust Stat Int
            std::vector<double> StempE(nExp, 0.0);
            subMatrix(RV_dot, VE_dot, nExp, nExp, nInt1, nExp, nInt1 + 1);
            matInv(VE_dot, nExp);
            for (size_t j = 0; j < nExp; j++) {
                for (size_t k = 0; k < nExp; k++) {
                    StempE[j] += (VE_dot[(nExp * j) + k] * beta_dot[k+1]);
                }
            }

            double statInt = 0.0;
            for (size_t j = 1; j < nExp1; j++) 
                statInt += beta_dot[j] * StempE[j-1];
            R_pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));

            // Robust Stat Joint
            std::vector<double> StempGE(nExp1, 0.0);
            subMatrix(RV_dot, A_i, nExp1, nExp1, nInt1, nExp1, 0);
            matInv(A_i, nExp1);
            for (size_t j = 0; j < nExp1; j++) {
                for (size_t k = 0; k < nExp1; k++) {
                    StempGE[j] += (A_i[(nExp1 * j) + k] * beta_dot[k]);
                }
            }

            double statJoint = 0.0;
            for (size_t k = 0; k < nExp1; k++) {
                statJoint += beta_dot[k] * StempGE[k];
	    }
            R_pvalJoint = (std::isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));
        }

        // Get model-based variance-covariance 
        for (size_t i = 0; i < nInt1 * nInt1; i++) {
            Vi_dot[i] *= sigma2;
        }

        // Model-based Interaction P-value
        std::vector<double> StempE(nExp, 0.0);
        subMatrix(Vi_dot, VE_dot, nExp, nExp, nInt1, nExp, nInt1 + 1);
        matInv(VE_dot, nExp);
        for (size_t j = 0; j < nExp; j++) {
            for (size_t k = 0; k < nExp; k++) {
                StempE[j] += (VE_dot[(nExp * j) + k] * beta_dot[k+1]);
            }
        }

        double statInt = 0.0;
        for (size_t j = 1; j < nExp1; j++) 
            statInt += beta_dot[j] * StempE[j-1];
        M_pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));


        // Model-based Joint P-value
        std::vector<double> StempGE(nExp1, 0.0);
        subMatrix(Vi_dot, A_i, nExp1, nExp1, nInt1, nExp1, 0);
        matInv(A_i, nExp1);
        for (size_t j = 0; j < nExp1; j++) {
            for (size_t k = 0; k < nExp1; k++) {
                StempGE[j] += (A_i[(nExp1 * j) + k] * beta_dot[k]);
            }
        }

        double statJoint = 0.0;
        for (size_t k = 0; k < nExp1; k++)
            statJoint += beta_dot[k] * StempGE[k];
        M_pvalJoint = (std::isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));


        if (subcategorical_exists) {
            for (size_t i = 0; i < nsubcategorical_columns; i++) {
                for (size_t j = 0; j < subcategorical_file_index[i].size(); j+=2) {
                    double n  = std::stod(values[subcategorical_file_index[i][j]]);
                    double af = std::stod(values[subcategorical_file_index[i][j+1]]);
                    strata_N[i]  += n;
                    strata_AF[i] += (n * ((!std::isnan(af)) ? af : 0));
                }
                strata_AF[i] = (strata_N[i] != 0) ? strata_AF[i] /= strata_N[i] : NAN;
            }
        }


        // Print variant info
        for (size_t ii = 0; ii < nVarInfo; ii++) {
            oss << values[varInfoIndices[ii]] << "\t";
        }

        if (subcategorical_exists) {
            for (size_t ii = 0; ii < nsubcategorical_columns; ii++) {
                oss << strata_N[ii] << "\t" << strata_AF[ii] << "\t";
            }
            std::fill(strata_N.begin(), strata_N.end(), 0.0);
            std::fill(strata_AF.begin(), strata_AF.end(), 0.0);			
        }

        oss << values[betaMarginalColumn] << "\t";
        if (robust) {
            oss << values[mbSeMarginalColumn] << "\t";
            oss << values[rbSeMarginalColumn ] << "\t";
        } else {
            oss << values[mbSeMarginalColumn] << "\t";
        }

        // Print betas
        for (size_t ii = printStart; ii < printEnd; ii++) {
            oss << beta_dot[ii] << "\t";
        }
        
        // Print robust variance-covariance first if applicable.
        if (robust) {
            for (size_t ii = printStart; ii < printEnd; ii++) {
                oss << sqrt(Vi_dot[ii * nInt1 + ii]) << "\t";
            }
            for (size_t ii = printStart; ii < printEnd; ii++) {
                for (size_t jj = printStart; jj < printEnd; jj++) {
                    if (ii < jj) {
                        oss << Vi_dot[ii * nInt1 + jj] << "\t";
                    }
                }
            }

            for (size_t ii = printStart; ii < printEnd; ii++) {
                oss << sqrt(RV_dot[ii * nInt1 + ii]) << "\t";
            }
            for (size_t ii = printStart; ii < printEnd; ii++) {
                for (size_t jj = printStart; jj < printEnd; jj++) {
                    if (ii < jj) {
                        oss << RV_dot[ii * nInt1 + jj] << "\t";
                    }
                }
            }
            
	    oss << values[mbPvalMarginalColumn ] << "\t" << M_pvalInt << "\t" << M_pvalJoint << "\t";
            oss << values[rbPvalMarginalColumn ] << "\t" <<  R_pvalInt << "\t" << R_pvalJoint << "\n";
        } else {
            for (size_t ii = printStart; ii < printEnd; ii++) {
                oss << sqrt(Vi_dot[ii * nInt1 + ii]) << "\t";
            }
            for (size_t ii = printStart; ii < printEnd; ii++) {
                for (size_t jj = printStart; jj < printEnd; jj++) {
                    if (ii < jj) {
                        oss << Vi_dot[ii * nInt1 + jj] << "\t";
                    }
                }
            }
            oss << values[mbPvalMarginalColumn] << "\t" << M_pvalInt << "\t" << M_pvalJoint << "\n";
        }
        p++;

        if (p % 100000 == 0) {
            results << oss.str();
            oss.str(std::string());
            oss.clear();
        }	
    }

    results << oss.str();
    oss.str(std::string());
    oss.clear();
    results.close();

    delete[] Vi_dot;
    delete[] VE_dot;
    delete[] O;
    delete[] O_dot;
    delete[] OV_prod;
    delete[] MRV_prod;
    delete[] RV_dot;
    delete[] A_i;

    cout << "Done.\n\nResults are in [" << cmd->outFile << "].";
    return;
}

void printTimeCompleted(double wall0, double wall1, double cpu0, double cpu1) 
{
    cout << "\n*********************************************************\n";
    cout << "Wall Time = " << wall1 - wall0 << " (sec)\n";
    cout << "CPU Time  = " << cpu1  - cpu0  << " (sec)\n\n";
}
