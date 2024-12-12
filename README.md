# REGEM

REGEM (RE-analysis of GEM summary statistics) is a software program for re-analyzing large-scale gene-environment interaction testing results, including multi-exposure interaction, joint, and marginal tests. It uses results directly from [GEM](https://github.com/large-scale-gxe-methods/GEM) output.

<br />  

Current version: 1.1

<br />  

## Contents 
- [Dependencies](#dependencies)
- [Quick Installation](#quick-installation)
- [Usage](#usage)
- [Contact](#contact)
- [License](#license)

<br />  

## Dependencies

- C++ compiler with C++11 support
- [Boost C++ Libraries](https://www.boost.org/) (Versions 1.70.0 - 1.79.0)
- Intel Math Kernal Library (MKL)

<br />

## Quick Installation

To install REGEM, run the following lines of code:
 ```
git clone https://github.com/large-scale-gxe-methods/REGEM
cd REGEM
cd src
make
 ``` 
 
<br />  

## Usage

### Running REGEM
- [Command Line Options](#command-line-options)
- [Input File](#input-file)
- [Output File](#output-file)
- [Examples](#examples)

<br />

### Command Line Options
Once REGEM is compiled, the executable ./REGEM can be used to run the program.
For a list of options, use ```./REGEM --help```.

<details>
     <summary> <b>List of Options</b> </summary>

```
General Options:
 
   --help                
   		Prints available options and exits.  
    
   --version             
   		Prints the version of REGEM and exits.


File Options:
   --input-file        
   		Path to the input file containing GEM results.  
    
   --out  
   		Full path and extension to where REGEM output results.  
   		Default: regem.out  
    
   --output-style  
   		Modifies the output of REGEM. Must be one of the following:
   			minimum: Output the summary statistics for only the GxE and marginal G terms.
   			   meta: 'minimum' output plus additional fields for the main G and any GxCovariate terms.
   				 For a robust analysis, additional columns for the model-based summary statistics will be included.
   			   full: 'meta' output plus additional fields needed for re-analyses of a subset of interactions.
   			Default: full  


Input File Options:  

   --exposure-names      
   		One or more column names in the input file naming the exposure(s) to be included in interaction tests.  
        
   --int-covar-names     
   		Any column names in the input file naming the covariate(s) for which interactions should be included for adjustment (mutually exclusive with --exposure-names).


Centering Conversion Options:

   --center-in
      A value of 0, 1, or 2 representing the original centering status in the input file: 0 indicates no centering, 1 indicates centering of all exposures and covariates, and 2 indicates centering of interaction covariates only.


   --center-out
      A value of 0, 1, or 2 representing the converted centering status: 0 for no centering, 1 to center all exposures and covariates, and 2 to center only the interaction covariates.


   --mean-value
      Pairs of variable names and their mean values, with each pair connected by an equal sign.

```
</details>

<br /> 

### Input File

REGEM takes as input an output file from a GEM run with the --output-style flag set to "full". The --output-style flag is available in GEMv1.4.1 and later versions of the software.

<br />

### Output File

REGEM will write results to the output file specified with the --out parameter (or 'regem.out' if no output file is specified).
Below are details of the possible column headers in the output file.

```diff 
SNPID              - The SNP identifier as retrieved from the genotype file.
CHR                - The chromosome of the SNP.
POS                - The physical position of the SNP. 
Non_Effect_Allele  - The allele not counted in association testing.  
Effect_Allele      - The allele that is counted in association testing.  
N_Samples          - The number of samples without missing genotypes.
AF                 - The allele frequency of the effect allele.
N_catE_*           - The number of non-missing samples in each combination of strata for all of the categorical exposures and interaction covariates.
AF_catE_*          - The allele frequency of the effect allele for each combination of strata for all of the catgorical exposure or interaction covariate.

Beta_Marginal           - The coefficient estimate for the marginal genetic effect (i.e., from a model with no interaction terms).
SE_Beta_Marginal        - The model-based SE associated with the marginal genetic effect estimate.  
robust_SE_Beta_Marginal - The robust SE associated with the marginal genetic effect estimate.

Beta_G             - The coefficient estimate for the genetic main effect (G).
Beta_G-*           - The coefficient estimate for the interaction or interaction covariate terms.
SE_Beta_G          - Model-based SE associated with the the genetic main effect (G).  
SE_Beta_G-*        - Model-based SE associated with any GxE or interaction covariate terms.
Cov_Beta_G_G-*     - Model-based covariance between the genetic main effect (G) and any GxE or interaction covariate terms.  
Cov_Beta_G-*_G-*   - Model-based covariance between any GxE or interaction covariate terms.
robust_SE_Beta_G   - Robust SE associated with the the genetic main effect (G).  
robust_SE_Beta_G-* - Robust SE associated with any GxE or interaction covariate terms.
robust_Cov_Beta_G_G-*   - Robust covariance between the genetic main effect (G) and any GxE or interaction covariate terms.
robust_Cov_Beta_G-*_G-* - Robust covariance between any GxE or interaction covariate terms. 

P_Value_Marginal           - Marginal genetic effect p-value from model-based SE.
P_Value_Interaction        - Interaction effect p-value (K degrees of freedom test of interaction effect) from model-based SE. (K is number of major exposures)
P_Value_Joint              - Joint test p-value (K+1 degrees of freedom test of genetic and interaction effect) from model-based SE.
robust_P_Value_Marginal    - Marginal genetic effect p-value from robust SE.
robust_P_Value_Interaction - Interaction effect p-value from robust SE.
robust_P_Value_Joint       - Joint test p-value (K+1 degrees of freedom test of genetic and interaction effect) from robust SE.
```

<br />

### Examples  

```unix
./REGEM --input-file gem.out --exposure-names cov1 --out regem_cov1.out
```

To change the centering status of cov1 from 0 to 1:

```unix
./REGEM --input-file gem.out --exposure-names cov1 --out regem_cov1.out --center-in 0 --center-out 1 --mean-value cov1=0.5
```
<br />
<br />

## Contact 
For comments, suggestions, bug reports and questions, please contact Han Chen (Han.Chen.2@uth.tmc.edu), Alisa Manning (AKMANNING@mgh.harvard.edu), or Kenny Westerman (KEWESTERMAN@mgh.harvard.edu). For bug reports, please include an example to reproduce the problem without having to access your confidential data.

<br />
<br />

## References
If you use REGEM, please cite
* Pham DT, Westerman KE, Pan C, Chen L, Srinivasan S, Isganaitis E, Vajravelu ME, Bacha F, Chernausek S, Gubitosi-Klug R, Divers J, Pihoker C, Marcovina SM, Manning AK, Chen H. (2023) Re-analysis and meta-analysis of summary statistics from gene-environment interaction studies. Bioinformatics 39(12):btad730. PubMed PMID: [**38039147**](https://www.ncbi.nlm.nih.gov/pubmed/38039147). PMCID: [**PMC10724851**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10724851/). DOI: [**10.1093/bioinformatics/btad730**](https://doi.org/10.1093/bioinformatics/btad730). 
 

<br />
<br />

## License 

 ```
 REGEM: RE-analysis of GEM summary statistics
 Copyright (C) 2021-2024 Duy T. Pham and Han Chen
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ```
