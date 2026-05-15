# MSMICA
**Mass Spectrometry Metabolomics Identification through Connections of Untargeted Metabolomics Data**

## Code Availability
___
**MSMICA** is available for non-commercial use only. All commercial use requires a business license. Please contact James Zhan at jamesjiadazhan@gmail.com for more information.

## License
___
MSMICA is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0). See the `LICENSE` file for details.

## Overview
___

**MSMICA** is an R package that uses experimentally optimized parameters — including adduct formations, adduct correlations, enzyme-based precursor-product correlations, and isotope correlations — to quantify the confidence of putative metabolite identification using a rigorous Bayesian probability framework.

The algorithm assigns each LC-MS feature a posterior identification probability and outputs a **unique** list of identified metabolites with one adduct form, m/z, and retention time per entry.

## Citation
___
If you use MSMICA in your research, **please cite our original work**:

- *Citation TBD — manuscript under review*

## How MSMICA Works
___

MSMICA performs identification in three main stages:

1. **m/z annotation** — Each observed feature is matched against the KEGG/HMDB metabolite database (combined or KEGG-only) at the user-specified ppm tolerance. All user-specified adduct forms are considered simultaneously.

2. **Clustering of adducts and isotopes** — For metabolites whose candidate features form adduct/isotopologue clusters (i.e. two or more adduct forms of the same metabolite co-elute and co-vary in intensity), a Bayesian model incorporating adduct formation, adduct correlation (Spearman r), and isotopologue abundance ratio is used to reduce the identity redundancies of features.  

3. **Local optimization per monoisotopic mass** — For all remaining features, MSMICA uses a local Bayesian optimization scores each candidate feature using retention-time prediction, precursor-product/transporter Spearman correlation, and biospecimen-specific HMDB concentration priors. This stage iterates until no new metabolite identification can be made.

## What MSMICA Can Do
___

**Input:**
- A LC-MS feature table with m/z, retention time, and per-sample intensities as columns (wide format).
- Optional: a list of pre-identified metabolites (e.g. from targeted analysis or authentic standards) with KEGG IDs, m/z, and retention times to serve as network anchors.
- Optional: a sample class file to restrict correlation analysis to study samples only.

**Output:**
- A unique list of identified metabolites with assigned adduct, m/z, retention time, identification method, and posterior probability (0–100 %).
- When applied to human plasma, urine, tissue, and mouse serum with combined HILIC positive / C18 negative data, MSMICA typically identifies ~3,000 unique metabolites per sample type.

## Dependency Installation
___

MSMICA is not yet on CRAN. Install from GitHub using **remotes** after installing the package dependencies below.

**Package dependencies:** `dplyr`, `readr`, `tidyr`, `data.table`, `mgcv`, `pracma`, `preprocessCore`, `imputeLCMD`, `MetaboCoreUtilsAdduct`.

**Installation helpers:** `remotes` for GitHub packages and `BiocManager` for Bioconductor packages.

```r
# Install only missing R package dependencies from CRAN
required_pkgs <- c("dplyr", "readr", "tidyr",
                   "data.table", "mgcv", "pracma", "remotes")
missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("preprocessCore", "pcaMethods", "impute"))

# Install imputeLCMD from CRAN
if (!requireNamespace("imputeLCMD", quietly = TRUE))
    install.packages("imputeLCMD")
```

## MSMICA Installation
```r
# Install MSMICA from GitHub
remotes::install_github("jamesjiadazhan/MSMICA")
```

To install from a local copy, install the dependencies above first. Then download and unzip the MSMICA repository, copy the local path, replace `"/path/to/MSMICA"` below, and run:

```r
install.packages("/path/to/MSMICA", repos = NULL, type = "source")
```

## Getting started
### Check the "Get started" header on the top of the MSMICA website
- ### [Click here](https://jamesjiadazhan.github.io/MSMICA_manual/articles/MSMICA.html)


## Contact
___
- **Author:** Jiada (James) Zhan — jzha832@emory.edu or jamesjiadazhan@gmail.com
- **Lab:** Dean Jones and Young-Mi Go's lab, Emory University, Atlanta, GA, USA
