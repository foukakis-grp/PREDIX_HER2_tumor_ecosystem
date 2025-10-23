# Proteomics Analysis (PREDIX Baseline)

## Author: Luca Gaessler
### Karolinska Institutet & SciLifeLab

Proteomics data was generated from patient biopsies of the clinical PREDIX trial using unlabeled data-independent acquisition (DIA) mass spectrometry. For details on data generation and biological interpretation, please refer to the main manuscript.

## Repository Structure
```
Proteomics/
├── data/              # Raw data files
├── figures/           # Output figures
├── processed/         # Processed data files
├── renv/              # R package management
├── scripts/           # Analysis scripts
├── Proteomics.Rproj   # R project
└── renv.lock          # R package versions
```

## R Environment
This project uses `renv` for reproducible R package management.

1. Open `Proteomics.Rproj` in RStudio
2. Install required packages:
```R
renv::restore()
```

## Contact
For questions about the proteomics analysis: luca.gaessler@ki.se