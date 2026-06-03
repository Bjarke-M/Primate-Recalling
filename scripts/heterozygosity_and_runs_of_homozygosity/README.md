# heterozygosity_and_runs_of_homozygosity

Per-sample heterozygosity per chromosome and PLINK-based detection of runs of homozygosity. Outputs CSV tables consumed by Figure 2.

## Files

- `workflow.py` — GWF workflow.
- `templates.py` — template functions used by the workflow.
- `het_stats.py` — heterozygosity calculation per base pair per chromosome.
- `roh_stats.py` — runs-of-homozygosity statistics from PLINK output.
- `primate_stats.yaml` — conda environment specification.

## Use

1. Create the environment: `conda env create -f primate_stats.yaml` and activate it.
2. From this directory, run `gwf run`.
