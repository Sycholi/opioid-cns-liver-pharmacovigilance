# Opioid-associated CNS depression in liver disease — analysis code

Analysis code for the manuscript *"Opioid-associated central nervous system depression in liver disease: triangulated evidence from two regulatory databases and a surgical cohort to guide nursing monitoring."*

The study triangulates two spontaneous reporting systems (FAERS, JADER) with clinical benchmarking in a real-world surgical cohort (INSPIRE) to quantify opioid-associated CNS-depression signals in liver disease: magnitude, agent ranking, time-to-onset, seriousness, and incidence/timing/severity in ward patients.

## Data sources

| Source | Access | Notes |
|---|---|---|
| FAERS/AERS quarterly ASCII files, 2004 Q1–2026 Q2 | Free download, https://fda.gov | Quarterly zip archives (~90 files) |
| JADER case reports, PMDA 202608 release | Free, https://www.pmda.go.jp | CP932-encoded CSVs inside a single zip |
| INSPIRE v1.4.2 | PhysioNet credentialed access + signed DUA | https://physionet.org/content/inspire — **not** redistributable; run on your own licensed copy |

No raw data are included in this repository.

## Requirements

- R 4.4.2 with `data.table`
- Python 3 with `matplotlib`, `pandas`
- A SLURM cluster is assumed by the submission scripts (`*.sh`); each `*.R` also runs standalone with `Rscript`

## Pipeline order

1. `gate_faers.sh` / `gate_faers.R` — download and parse FAERS quarters; build report-level tables (dedup by case version, primary-suspect roles)
2. `faers_analysis.R` — overall and per-drug disproportionality (ROR with Woolf CI, PRR with Evans criteria, IC025 beta-binomial), Weibull time-to-onset, liver stratum, subgroup stratification
3. `sensitivity_analysis.R` + `sens2_analysis.R` — sensitivity suite: liver-definition variants (primary union / indication-only / history-term-only), expanded nine-term CNS set, secondary-suspect roles, report-level outcome proportions, JADER t≥0 rerun
4. `jader_analysis.R` / `jader_gate.R` — JADER aligned five-opioid set: ROR, outcomes, Weibull TTO
5. `inspire_validation.R` — INSPIRE cohort: postoperative opioid recipients with baseline ward GCS ≥13, incidence of GCS ≤12, severity, timing, MELD-Na stratification (outputs aggregate statistics only)
6. `make_figures.py` / `make_figures4.py` — all main-text figures (forest plots, Weibull TTO panels, outcome proportions, INSPIRE benchmarking)

Paths assume a working directory `~/nurse_faers/` holding the FAERS/JADER downloads with outputs in `~/nurse_faers/analysis/`; adjust to your environment.

## Outputs ↔ manuscript mapping

| Script output | Manuscript element |
|---|---|
| `faers_*ror.csv`, `sens2` liver CSVs | Tables 2–3, Figures 1B, 2 |
| `*_tto.csv`, Weibull fits | Figure 3A |
| `*_outcome*_reports.csv` | Figure 3B |
| `inspire_*` aggregate tables | Table 4, Figure 4 |

## Citation

If you use this code, please cite the manuscript (Li Y, Liu Y, et al., *International Journal of Nursing Sciences*, in submission 2026).

## Authors

- Yu Li — Department of Hepatobiliary Surgery, Sun Yat-sen University Cancer Center, Guangzhou, China
- Yuyao Liu (corresponding) — Department of Oncology, The First Affiliated Hospital of Jinan University, Guangzhou, China

## License

MIT
