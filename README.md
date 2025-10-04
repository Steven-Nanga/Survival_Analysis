# Survival_Analysis
Collection of survival analyis mini projects
## Survival analysis scripts

This project uses classic survival datasets from the R `survival` package (see below for details). All data files are in the `data/` folder. Analysis scripts are in `scripts/`, and outputs are saved in `outputs/`.

To run an analysis, use:
```powershell
Rscript scripts/analyze_colon.R
Rscript scripts/analyze_lung.R
Rscript scripts/analyze_ovarian.R
Rscript scripts/analyze_stanford2.R
Rscript scripts/analyze_veteran.R
Rscript scripts/analyze_kidney.R
```

### Dataset explanations
- `colon.csv`: Colon cancer trial data (`colon` from R survival package)
- `lung.csv`: Advanced lung cancer data (`lung` from R survival package)
- `ovarian.csv`: Ovarian cancer data (`ovarian` from R survival package)
- `stanford2.csv`: Heart transplant data (`stanford2` from R survival package)
- `veteran.csv`: Veteran's lung cancer trial (`veteran` from R survival package)
- `kidney.xlsx`: Kidney infection recurrence data (`kidney` from R survival package)

For more details, see `?colon`, `?lung`, etc. in R after loading the `survival` package.
