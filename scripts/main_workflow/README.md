# main_workflow

GWF workflow covering the full pipeline: FASTQ download, paired-end short-read mapping, per-individual variant calling, and joint genotyping.

## Files

- `workflow.py` — GWF target definitions, organised in blocks.
- `templates.py` — template functions invoked by `workflow.py`.
- `env.yaml` — conda environment specification.

## Use

1. Create the environment: `conda env create -f env.yaml` and activate it.
2. From this directory, run `gwf run`.

Run a few blocks at a time rather than the whole workflow at once.
