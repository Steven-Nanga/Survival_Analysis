Survival analysis scripts for classic datasets from the R `survival` package

## Folder structure
- `data/`: Datasets (from R survival package)
- `scripts/`: R analysis scripts
- `outputs/`: Results and plots

## What each analysis script does
- Loads the dataset
- Cleans and summarizes demographics
- Plots Kaplan-Meier curves (with log-rank test and risk tables)
- Fits univariable and multivariable Cox proportional hazards models
- Checks proportional hazards assumptions (Schoenfeld residuals)
- Exports publication-quality tables (gtsummary, Word)
- Saves outputs to `outputs/<dataset>/`

## Requirements
- R >= 4.0
- Packages: survival, survminer, ggplot2, dplyr, tidyr, gtsummary, broom, flextable, officer, MASS, readr, readxl

## How to run
From the repository root in PowerShell run:
```powershell
Rscript scripts/analyze_colon.R
Rscript scripts/analyze_lung.R
Rscript scripts/analyze_ovarian.R
Rscript scripts/analyze_stanford2.R
Rscript scripts/analyze_veteran.R
Rscript scripts/analyze_kidney.R
```

## Dataset explanations
- `colon.csv`: Colon cancer trial data (`colon` from R survival package)
- `lung.csv`: Advanced lung cancer data (`lung` from R survival package)
- `ovarian.csv`: Ovarian cancer data (`ovarian` from R survival package)
- `stanford2.csv`: Heart transplant data (`stanford2` from R survival package)
- `veteran.csv`: Veteran's lung cancer trial (`veteran` from R survival package)
- `kidney.xlsx`: Kidney infection recurrence data (`kidney` from R survival package)

For more details, see `?colon`, `?lung`, etc. in R after loading the `survival` package.
