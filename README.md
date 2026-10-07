# Opioid-associated CNS depression in liver disease — analysis code

Analysis code for the manuscript *"Opioid-associated central nervous system depression in liver disease: triangulated evidence from two regulatory databases and a surgical cohort to guide nursing monitoring."*

The study triangulates FAERS and JADER spontaneous reporting data with clinical benchmarking in the INSPIRE surgical cohort to quantify opioid-associated CNS-depression signals in liver disease.

## Contents

`analysis/` contains the R and Python scripts for data parsing, disproportionality analysis (ROR/PRR/IC025), Weibull time-to-onset modelling, sensitivity analyses, the INSPIRE cohort analysis, and all main-text figures. Requirements: R 4.4.2 (`data.table`) and Python 3 (`matplotlib`, `pandas`).

## Data

No raw data are included. FAERS (fda.gov) and JADER (pmda.go.jp) are freely downloadable; INSPIRE requires PhysioNet credentialed access and a signed data use agreement (no redistribution).

## License

MIT
