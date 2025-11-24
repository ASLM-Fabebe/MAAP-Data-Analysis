# Introduction

These scripts have been developed to simplify and standardize AMR/C/U data analysis through ASLM as part of the MAAP Phase II project funded through Fleming Fund. These scripts are largely based on the free and open-source [AMR package](https://msberends.github.io/AMR/index.html) implemented in R.

# Set A: AMR analysis script- MAAP2 (ASLM)
Follow the steps below:

## Step 1: Set-Up and Software Installation
To get started, install the required software:
- Install R: Visit [cran.r-project.org](https://cran.r-project.org/bin/windows/base/) and download the appropriate version for your system (e.g., Download R-4.4.3 for Windows or Mac).
- Install RStudio: Go to the [RStudio desktop download page](https://posit.co/download/rstudio-desktop/) and download the latest version. Click the installer to begin installation once downloaded.
- Install GitHub Desktop: Navigate to the [GitHub-Desktop download page](https://desktop.github.com/download/) and click Download for Windows or Mac. Follow the installation instructions and open the application after successful installation.

## Step 2: Clone the Repository
Clone the repository to your local machine using GitHub Desktop. The repository contains scripts for AMR analysis, including AMC & AMU processing and plotting.
- Go to `File` -\> `Clone repository` -\> URL and then paste the following link `https://github.com/ASLM-Fabebe/MAAP-Data-Analysis.git` in the first box and
- Click Clone
- In RStudio go to `File` -\> `New Project` -\> `Existing Project` -\> Browse (Navigate to your Documents -\> GitHub -\> MAAP-Data-Analysis)
- Click MAAP-Data-Analysis.Rproj
- This will open the MAAP-Data-Analysis scripts within RStudio on your local computer.

## Step 3: Upload your Dataset
- Add your input data to the `test-data` folder. For now the script only accepts input in Excel format. Importantly, the file with AST data should have "AMR" as its prefix.
- Now you can attempt to run the analysis script.

## Step 4: Running the Analysis
Execute the main script below - this will open an app that will walk you through the data preparation process, executing a Play-along mode.

```{r}
library(shiny); runApp('amr_analysis_dev.R')
```

- You will then need to enter your country name and register it by clicking the 'Register Country' button.
- Then you will need to match your datasets' variables with those of the system for alignment. Once this step is done, you will need to save by clicking the 'Save Data' button. The file will be saved in MAAP-Data-Analysis/test-data/analysis_update.

Finally, you will click the 'Begin Analysis' button. This will perform the analysis process end-to-end automatically.

## Output
If everything runs successfully, you should have the following in the results folder:

1. 4 tables based on the provided input file
- Demographics
- Facilities information
- Generic organisms list
- The test result file (interpreted AST results), with cleaned and standardized AST interpretations
  
2. An antibiogram of the tested bug-drug combinations

3. Sub-folders of analyzed pathogens e.g. ESKAPE pathogens, with each folder named after the organism analyzed, the results in these folders include:
- CSV files of Resistance distribution by Age, gender, and specimen type
- Barplots of the resistance prevalence of ESKAPE pathogens
 
# Set B: AMC analysis script- MAAP2 (ASLM)
This R script processes and analyzes antimicrobial consumption (AMC) data to calculate Defined Daily Doses (DDD) and DDD per 1,000 inhabitants per day (DiD), with visualization of trends and patterns. Please follow the Steps 1 to 4 below to be able to utilize the scripts develop.

## Step 1: Set-Up and Software Installation
To get started, install the required software:
- Install R: Visit [cran.r-project.org](https://cran.r-project.org/bin/windows/base/) and download the appropriate version for your system (e.g., Download R-4.4.3 for Windows or Mac).
- Install RStudio: Go to the [RStudio desktop download page](https://posit.co/download/rstudio-desktop/) and download the latest version. Click the installer to begin installation once downloaded.
- Install GitHub Desktop: Navigate to the [GitHub-Desktop download page](https://desktop.github.com/download/) and click Download for Windows or Mac. Follow the installation instructions and open the application after successful installation.

## Step 2: Clone the Repository
If you haven't already, clone the repository to your local machine using GitHub Desktop. Follow the steps up in the AMR module (Set A). The repository contains scripts for AMR analysis, including AMC & AMU processing and plotting.

## Step 3: Upload Input Data
Place your input datasets into the test-data/AMC/ folder. Two datasets are required:

1. Reference Data
- Format: .xlsx (Excel)
- File name: ATC-DDD WHO core and optional antimicrobials.xlsx
- Description: Core and optional antimicrobial reference list from WHO.

2. AMC Test Data
- Format: .xlsx (Excel)
- File name: AMC_test_data.xlsx (only for demo)
- Description: Country-specific antimicrobial consumption data, including product names, strengths, routes, and quantities.

### Importantly, the dataset file should have "AMC" as its prefix.

## Step 4: Run the Analysis
To run the analysis in RStudio:
- Open RStudio and go to: File → New File → R Script
- Copy and paste the single line of command below into the new script:

```{r}
library(shiny); runApp('amc_analysis_app.R')
```
This will launch an app that will walk you through a 6 step process for data preparation, cleaning, and finally analysis.

- The first step, you’ll be prompted to enter your country’s name and estimated population.
-
- The script will then execute an end-to-end analysis automatically.

## Output
If the script runs successfully, the following visualizations will be generated and saved in the plots_AMC folder:
- AMC Class Distribution (percentage)
- AMC Class Totals
- AWaRe Category Distribution
- Monthly Consumption Trends

## Usage
1. Set the population size when prompted
2. Ensure input files are in the correct folder (`test-data/AMC/`)
3. Run the script to generate:
   - Processed data frame (`amc`) with DiD calculations
   - Visualizations in `plots_AMC/` directory

# Set C: AMU analysis script- MAAP2 (ASLM)
This script analyzes Antimicrobial Use (AMU) data from healthcare facilities, producing:
- Demographic summaries (age/gender distributions)
- AMU prevalence metrics (overall, by ward/age/gender)
- Antimicrobial utilization patterns (by class/molecule/route/AWaRe category)

```{r}
library(shiny); runApp('amu_app_dev.R')
```
## Outputs Created
If the script runs successfully, the following outputs will be generated and saved in the plots_AMU folder:
1. Tables (Excel Files)
- Patient % by age group per facility
- Patient % by gender per facility
- Overall AMU prevalence % per facility
- AMU % by age group per facility
- AMU % by gender per facility
- Avg. antimicrobials per patient per facility
- AMU % by ATC class per facility
- AMU % by antimicrobial molecule per facility
- AMU % by AWaRe category per facility
- AMU % by administration route per facility
- 
2. Plots (PNG Files)
- AMU Prevalence by Site
- AMU by ATC Class
- AMU by AWaRe Category
- AMU by Route


