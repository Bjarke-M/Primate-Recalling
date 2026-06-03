# additional_scripts

Helper scripts used by the main workflow.

## Scripts

- `download_fastq_per_individual.py` — download FASTQ files for one individual.
- `download_fastq_per_individual_check_md5.py` — same, with md5 verification.
- `check_md5.py` — verify md5 checksums of downloaded files.
- `make_reference.py` — prepare a reference genome (index, dictionary).
- `make_regions_file.py` — build per-reference region/interval files for parallel calling.
- `find_chrX.py` — locate X-chromosome contigs in a reference.
- `make_geno_metadata.py` — generate genotyping metadata.
- `make_loc_metadata.py` — generate sample-location metadata.
- `make_genDB_folder_and_map.py`, `make_genDB_folder_and_map2.py` — set up GenomicsDB folders and sample maps for joint genotyping.
- `make_simplified_batch_file.py` — simplify batch input files for GWF targets.

## Files

- `special_contigs/` — manually curated contig lists for non-standard references.
