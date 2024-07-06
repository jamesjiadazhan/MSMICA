# MSMICA
**Mega scale putative metabolomics identification through connections of untargeted metabolomics data** 

## Code availability
___
**MSMICA** is available to only non-commercial use. All commercial use of MSMICA requires a business license. Please contact Jiada (James) Zhan at jzha832@emory.edu or the Office of Technology Transfer at Emory University at mcoburn@emory.edu for more information.

## Overview
___

**MSMICA** is an R package that uses experimentally optimized parameters, including adduct formations, adduct correlations, enzyme-based precursor-product correlations, and isotope correlations, to quantify the confidence of putative metabolite identification with a rigorous Bayesian probability approach. 


## User-friendly tutorial page
___
- ### [Click here]()


## Citation
___
If you are using the MSMICA in your research, **please be sure to cite our original work**. By doing so, you not only add credibility to your findings but also recognize and appreciate our intellectual efforts and contributions. The appropriate citation is as follows:

- *bioRxiv citation or Nature Methods citation*
  - website link

## How MSMICA works
___
The **MSMICA** package performs calculations in 3 steps:

1. m/z annotation using KEGG database using users' provided adduct forms (default 5 ppm).

2. Bayesian probability calculation for metabolites with clustering patterns using adduct correlations, isotope correlations and adduct formations.

3. Bayesian probability calculation for all the other metabolomics features using precursor-product correlations, using KEGG precursor-product enzyme-based reaction database.

## What MSMICA can do
___

- Input: 
    - The algorithm starts with the input of a liquid chromatography-mass spectrometry (LC-MS) feature table with a list of features’ m/z, retention time, and intensities. 
        - The default is an untargeted approach. 
        - Users have the option of providing a list of identified or confidently annotated metabolites as additional algorithm starting points with column mode, metabolite names, KEGGID, m/z, retention times, and intensities. This may further increase the coverage of MSMICA.

- Output: 
    - The output is a unique list of identified metabolites with the one adduct form, m/z, and retention time assigned. 
        - Regarding output, by applying MSMICA to studies with human plasma, human urine, human tissue, and mouse serum samples, we identified about 3000 unique metabolites by combining HILIC positive and C18 negative results in each sample type. 
        - Using our retention time reference library (m/z 5 ppm, RT 30 seconds), we confirmed some of these metabolites with an average fraction correct rate of 0.96.

## Installation
___

Currently, **dietaryindex** is not available on [CRAN]

To install MSMICA from GitHub, use the **devtools** package:

Package dependencies: **dplyr**, **readr**, **tidyr**, **progress**, **MetaboCoreUtilsAdduct**.

```
# If you don't have the following dependencies installed already
install.packages("dplyr")
install.packages("readr")
install.packages("tidyr")
install.packages("progress")
install.packages("devtools")
devtools::install_github("jamesjiadazhan/MetaboCoreUtilsAdduct") # Install the package from GitHub

# Now, install MSMICA
devtools::install_github("jamesjiadazhan/MSMICA") # Install the package from GitHub
```

To install MSMICA locally, you can first download the package manually by clicking the **Code** button and downloading the ZIP file. Unzip the file. Then, use the following codes (replace my path with your path):

```
# If you don't have the following dependencies installed already
install.packages("dplyr")
install.packages("readr")
install.packages("tidyr")
install.packages("progress")
install.packages("devtools")
devtools::install_github("jamesjiadazhan/MetaboCoreUtilsAdduct") # Install the package from GitHub

# Now, install MSMICA
install.packages("/Users/yan/Downloads/MSMICA-main", repos = NULL, type = "source")
```


If something happens like the following, first try to enter 1 in the terminal (lower box). If not successful, then try to enter 2. **It will take a while if you are a new R user.**
```
  These packages have more recent versions available.
  It is recommended to update all of them.
  Which would you like to update?

  1: All                          
  2: CRAN packages only           
  3: None                         
  4: tzdb  (0.3.0 -> 0.4.0) [CRAN]
  5: vroom (1.6.1 -> 1.6.3) [CRAN]
```


## Getting Started
___
- ### [Click here for codes with printing outputs]()
```
# Load the necessary dependency packages for MSMICA
## this is for data processing
library(dplyr)
## this is for tibble data structure
library(readr)
## this is for data processing
library(tidyr)
## this is for progress bar display
library(progress)
## this is for calculating theoretical m/z based on the KEGG monoisotopic masses and adduct forms
library(MetaboCoreUtilsAdduct)

# Load MSMICA
library(MSMICA)

# The example data includes a feature table from a LC-MS metabolomics study involving 7 paired kidney cancer tissue samples and normal kidney tissue samples taken from the same subjects (7 subjects, 14 actual samples, 42 triplicate samples). These samples were analyzed using HILIC positive. The feature table is a data frame with the first column as m/z (mass-to-charge ratio), the second column as retention time, and the rest of the columns as the intensity values for each sample.
data(feature_table_exp_hilicpos)

# Set a working directory for the MSMICA output files. This is where the output files will be saved.
setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PhD dissertation/MSMICA/Publication/Abstract/MSMICA algorithm/MSMICA example data")

# Filter out features appear less than 20% of all samples. The intensity column starts from the 3rd column.
feature_table_exp_hilicpos = QC_filter(x = feature_table_exp_hilicpos, metabolite_start_column = 3, minimum_sample_appear = 0.20)

# Run the MSMICA algorithm with the example data. Remember to select only one ion mode at a time and select the appropriate adduct forms.
MSMICA_algorithm(met_raw_wide = feature_table_exp_hilicpos, Adduct = c("M+H","M+2Na-H","M+Na","M-H2O+H","M+K","M+2H"), prefix="MSMICA_test_hilicpos", ion_mode="positive")

```