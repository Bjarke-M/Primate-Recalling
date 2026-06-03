# Primate-Recalling

Supplementary code for:

**A Curated Genome-Scale Nucleotide Diversity Panel of Non-Human Primates**
Vasili Pankratov¹‡, Bjarke Meyer Pedersen¹‡, Erik Fogh Sørensen¹, Kasper Munch¹, Thomas Bataillon¹, Mikkel Heide Schierup¹, Juraj Bergman¹\*

¹ Bioinformatics Research Centre, Aarhus University, DK-8000 Aarhus C, Denmark
‡ Equal contribution
\* Corresponding author: jurajbergman@birc.au.dk

## Repository layout

- `scripts/` — workflows and helper scripts for downloading reads, mapping, variant calling, joint genotyping, sex-chromosome annotation, and heterozygosity/ROH analysis.
- `metadata/` — supplementary tables: sample statistics, sequencing accessions, and reference genome sources.
- `plots/` — notebooks, R scripts, and input data used to produce the figures in the paper.

Each subdirectory contains its own README with details on contents and use.

## Reproducing the pipeline

1. Set up the environment in `scripts/main_workflow/env.yaml`.
2. Run the main GWF workflow in `scripts/main_workflow/` to produce BAMs and joint-genotyped VCFs.
3. Run the per-sample analyses in `scripts/heterozygosity_and_runs_of_homozygosity/` and `scripts/annotate_sex_chr/`.
4. Use the notebooks in `plots/` to regenerate the figures from the resulting tables.
