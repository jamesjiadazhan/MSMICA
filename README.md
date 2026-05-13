# MSMICA
**Mega-Scale Metabolomics Identification through Connections of Untargeted Metabolomics Data**

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

## Installation
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

# Install MetaboCoreUtilsAdduct from GitHub
remotes::install_github("jamesjiadazhan/MetaboCoreUtilsAdduct")

# Install MSMICA from GitHub
remotes::install_github("jamesjiadazhan/MSMICA")
```

To install from a local copy, install the dependencies above first. Then download and unzip the MSMICA repository, copy the local path, replace `"/path/to/MSMICA"` below, and run:

```r
install.packages("/path/to/MSMICA", repos = NULL, type = "source")
```

If you receive MSMICA as a local source archive, such as `.tar.gz`, install dependencies first and then use:

```r
install.packages("/path/to/MSMICA_1.0.0.tar.gz", repos = NULL, type = "source")
```

## Getting Started
___

```r
library(MSMICA)

# ── Step 1: Load your feature table ──────────────────────────────────────────
# The feature table must have m/z as column 1, retention time (seconds) as column 2, and per-sample intensities in the remaining columns.
data(feature_table_exp_hilicpos)

# ── Step 2: QC filter ─────────────────────────────────────────────────────────
# Remove features that appear in fewer than 20 % of samples.
feature_table_exp_hilicpos <- QC_filter(
    x = feature_table_exp_hilicpos,
    metabolite_start_column = 3,
    minimum_sample_appear = 0.20
)

# ── Step 3: Run MSMICA ────────────────────────────────────────────────────────
# Select one ion mode at a time and provide the appropriate adduct list.
adducts <- msmica_adducts(mode = "positive", sample_type = "fluid")

MSMICA_algorithm(
    met_raw_wide    = feature_table_exp_hilicpos,
    LC              = "HILIC",    # chromatography
    LC_run_time     = 5,          # minutes
    mz_threshold    = 10,
    ion_mode        = "positive",
    All_Adduct      = adducts,
    biospecimen     = "Blood",
    reaction_database = c("mammalia"),
    prefix          = "MSMICA_test_hilicpos"
)
```


## Function Reference
___

| Function | Description |
|---|---|
| `MSMICA_algorithm()` | Main entry point. Runs the full three-stage identification pipeline. |
| `QC_filter()` | Removes low-prevalence features (appear in < x% of samples). |
| `msmica_adducts()` | Returns preset adduct vectors for common ion mode and sample type combinations. |
| `find.Overlapping.mzs()` | Fast ppm-based m/z matching between two feature tables using `data.table`. |
| `custom_biochemical_reaction_loading()` | Loads the bundled curated biochemical reaction dataset. |

All other functions in the package are internal helpers called automatically by `MSMICA_algorithm()`.

## Key Parameters
___

| Parameter | Default | Description | Alternative Options |
|---|---|---|---|
| `mz_threshold` | `10` | m/z matching tolerance in ppm. Use 10 ppm with high-resolution instruments (including Orbitrap MS). | User-defined numeric threshold based on instrument performance. |
| `LC` | `"HILIC"` | LC column type for RT prediction. | `"RP"` or `"C18"`. |
| `LC_run_time` | — | Total LC run time in **minutes** (required). | Any positive numeric runtime in minutes. |
| `biospecimen` | `"Blood"` | Biospecimen type used for HMDB concentration priors. Only Blood and Urine are well characterized. Other biospecimen types are not well documented so using Blood may be a good alternative. | `"Urine"`, `"Feces"`, `"Cerebrospinal Fluid"`, `"Saliva"`, `"Breast Milk"`, `"Sweat"`, `"Cellular Cytoplasm"`, `"Amniotic Fluid"`, `"Aqueous Humour"`, `"Ascites Fluid"`, `"Lymph"`, `"Tears"`, `"Bile"`, `"Semen"`, `"Pericardial Effusion"`. |
| `ion_mode` | `"positive"` | Ionization mode. | `"negative"`. |
| `All_Adduct` | `msmica_adducts("positive", "fluid")` | Adduct forms considered for matching. | Use `msmica_adducts(mode, sample_type)` for presets: `mode = "positive"` or `"negative"`; `sample_type = "fluid"` or `"tissue"`. You can also provide a custom character vector. |
| `adduct_correlation_r_threshold` | `0.39` | Spearman correlation threshold for adduct correlation analysis. | User-defined numeric threshold (typically between 0 and 1). |
| `adduct_correlation_time_threshold` | `6` | Retention-time threshold (seconds) for adduct correlation analysis. | User-defined positive numeric value in seconds. |
| `isotopic_correlation_r_threshold` | `0.71` | Spearman correlation threshold for isotopic correlation analysis. | User-defined numeric threshold (typically between 0 and 1). |
| `isotopic_correlation_time_threshold` | `4` | Retention-time threshold (seconds) for isotopic correlation analysis. | User-defined positive numeric value in seconds. |
| `reaction_database` | `"mammalia"` | Biochemical reaction database(s) for precursor-product scoring. | `"general"`. |
| `imputation_method` | `"half_min"` | Missing-value imputation method. | `"QRILC"` or `NA` (no imputation). |
| `detail` | `FALSE` | Save intermediate CSVs (warning: 10+ large files with hundreds of MB). | `TRUE`. |
| `progress_log` | `FALSE` | Write all messages to a `.txt` log file. | `TRUE`. |

## Adduct Presets
___

`msmica_adducts()` returns preset adduct vectors for common ion mode and sample type combinations. These presets are meant as convenient starting points; users can still provide a custom character vector to `All_Adduct`.

| Preset | Exact adduct vector |
|---|---|
| `msmica_adducts("positive", "fluid")` | `c("M+H", "M+Na", "M+2Na-H", "M+H-H2O", "M+H-NH3", "M+ACN+H", "M+ACN+2H", "2M+H", "M+2H", "M+H-2H2O")` |
| `msmica_adducts("negative", "fluid")` | `c("M-H", "M+Cl", "M+FA-H", "M+Hac-H", "M-H+HCOONa", "M+Na-2H", "M-2H", "2M-H", "M+ACN-H")` |
| `msmica_adducts("positive", "tissue")` | `c("M+H", "M+K", "M+2K-H", "M+H-H2O", "M+H-NH3", "M+ACN+H", "M+ACN+2H", "2M+H", "M+2H", "M+H-2H2O")` |
| `msmica_adducts("negative", "tissue")` | `c("M-H", "M+Cl", "M+FA-H", "M+Hac-H", "M-H+HCOOK", "M-2H", "2M-H", "M+ACN-H", "M+K-2H")` |



## Contact
___
- **Author:** Jiada (James) Zhan — jzha832@emory.edu or jamesjiadazhan@gmail.com
- **Lab:** Dean Jones and Young-Mi Go's lab, Emory University, Atlanta, GA, USA
