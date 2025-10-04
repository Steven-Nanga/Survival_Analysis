Survival analysis scripts for the `data/` datasets

What this folder contains
- `00_setup_utils.R`: helper functions for package installation and folder creation
- `analyze_<dataset>.R`: analysis scripts for each dataset (colon, lung, ovarian, stanford2, veteran, kidney)

What each analysis script does
- Loads the dataset
- Basic data cleaning and demographics table
- Kaplan-Meier curves (with log-rank test and ggplot2 output)
- Univariable and multivariable Cox proportional hazards models
- PH assumption checks (Schoenfeld residuals)
- Publication-quality plots saved to `outputs/<dataset>/`

Requirements
- R >= 4.0
- Packages: survival, survminer, ggplot2, dplyr, tidyr, gtsummary, broom, flextable, officer, MASS, readr, readxl

Run
From the repository root in PowerShell run:
```powershell
Rscript scripts/analyze_colon.R
Rscript scripts/analyze_lung.R
Rscript scripts/analyze_ovarian.R
Rscript scripts/analyze_stanford2.R
Rscript scripts/analyze_veteran.R
Rscript scripts/analyze_kidney.R
```

Outputs
- Plots and tables are written to `outputs/<dataset>/` as PNG and RDS files.

Notes
- These scripts are written for the CSV/XLSX formats present in `data/`. If your data columns differ, update variable names near the top of each script.
