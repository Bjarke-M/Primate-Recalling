# scripts

Pipelines and helper scripts.

- `main_workflow/` — GWF workflow: fastq download, mapping, per-sample variant calling, joint genotyping.
- `annotate_sex_chr/` — identify sex-chromosome-derived contigs in references and assign genetic sex (XX/XY) to each sample.
- `heterozygosity_and_runs_of_homozygosity/` — per-sample heterozygosity per chromosome and PLINK-based runs of homozygosity.
- `additional_scripts/` — utilities used by the workflows (reference preparation, region files, metadata generation, md5 checks).
