# main_workflow

GWF workflow covering the full pipeline: reference preparation, FASTQ download, paired-end short-read mapping, per-individual variant calling, and joint genotyping per species.

## Files

- `workflow.py` — GWF target definitions, organised in four blocks (A–D).
- `templates.py` — template functions invoked by `workflow.py`.
- `env.yaml` — conda environment specification.

## Use

1. Create the environment: `conda env create -f env.yaml` and activate it.
2. Place `supplementary_table_2.txt` and `supplementary_table_3.txt` (from `../../metadata/`) where the paths at the top of `workflow.py` expect them.
3. From this directory, run `gwf run`.

Run the blocks in order (A → B → C → D), and only a few blocks at a time. Each step writes a sentinel file under `<group>/done/` that the next step depends on; this lets the workflow resume after interruption and skip already-finished work.

## Block A — Reference preparation

For each reference genome listed in `supplementary_table_3.txt`:

- **A.0** Create the `ref/`, `fastq/`, `bam/`, `gVCF/`, and `done/` directories for the reference.
- **A.1 `download_ref`** Download the reference FASTA from its FTP URL.
- **A.2 `mask_reference`** *(optional, commented out)* Mask regions listed in `masked_regions.bed`. Enable by swapping the commented input paths in A.3–A.5.
- **A.3 `cut_contigs`** Drop contigs shorter than 1000 bp. 
- **A.4 `make_fasta`** Build the `.fai` index and sequence dictionary (`samtools faidx`, `picard CreateSequenceDictionary`) required by BWA and GATK.
- **A.5 `make_regions`** Generate a per-reference regions table. This table is manually edited after the run to set per-region ploidy for X, Y, PAR, and mitochondrial contigs.

## Block B — FASTQ download and uBAM preparation

For every individual in `supplementary_table_2.txt`:

- **B.1 `download_per_individual`** Download all FASTQ files for the individual in one job, with md5 verification against the values in the metadata table.
- **B.2 `concatfastqs` / `renamefastqs`** Individuals with multiple runs have their R1/R2 files concatenated; single-run individuals are renamed to the same `{ind}_R1.fastq.gz` / `{ind}_R2.fastq.gz` convention so downstream jobs need only one code path.
- **B.3 `makeuBAM`** Convert paired FASTQs to an unmapped BAM (`picard FastqToSam`), preserving read-group metadata for GATK.
- **B.4 `splituBAM`** Split the uBAM into shards. The number of shards is written to `{ind}_nsplitubams.txt`, which Block C reads to fan out mapping.

## Block C — Mapping

Mapping is sharded per individual; Block C only schedules jobs once `splituBAM` has finished and `{ind}_nsplitubams.txt` exists.

- **C.1 `markadapt`** Mark sequencing adapters in each shard (`picard MarkIlluminaAdapters`) so BWA soft-clips them.
- **C.2 `mapBAM`** Map each shard with BWA-MEM against the reference produced in Block A, merging mapping metadata back from the uBAM (`picard MergeBamAlignment`) in a single piped command.
- **C.3 `mergeBAMs`** Merge all mapped shards for an individual into one BAM.
- **C.4 `markduplicates`** Mark and remove PCR/optical duplicates (`picard MarkDuplicates`).
- **C.5 `coordsort`** Coordinate-sort and index the final BAM.
- **C.6 `cov`** Compute coverage per region using the regions table from A.5. Used downstream to assign genetic sex and to flag samples that fail coverage thresholds.
- **C.7 `cov_batched` + `concatenate_cov_files`** *(alternative, commented out)* For references with very many contigs, coverage is computed in batches of 10 000 regions and concatenated. Swap C.6 for C.7 when the single-job version exceeds wall-time.

## Block D — Calling and genotyping

Calling happens per individual; genotyping is joint per species (`GVCF_FOLDER`).

Pre-calling:

- **D.1 `find_chrX`** Identify X-chromosome contigs (≥ 1 Mb) in each reference and emit `regions_<ref>_updated.txt` with the X annotation merged in.
- **D.2 `make_simplified_batch_file`** Collapse the regions table into a `batches_<ref>.txt` file that groups regions into calling batches with their female/male ploidy.

Per-individual calling:

- **D.3 `make_batch_metadata`** Build per-individual BED files and GATK `.intervals` files for each (batch, ploidy) combination, using each sample's genetic sex from `samples_coverage_stats.txt`. Regions longer than 2 Mb are split into chunks.
- **D.4 `call_batch_with_bed`** Run `gatk HaplotypeCaller` in GVCF mode for every (individual, batch, ploidy) combination. Ploidy is set per region so that female chrX is called diploid, male chrX/Y haploid, mitochondrial DNA haploid, etc.

Joint genotyping per species:

- **D.5 `make_geno_metadata`** Build the genotyping interval files for the cohort, again chunked at 2 Mb.
- **D.6 `make_genDB_folder_and_map`** Create a GenomicsDB workspace folder and `cohort.sample_map` for each (batch, female_ploidy, male_ploidy) combination. Three cases are handled explicitly:
  - matched ploidies (autosomes, PAR, mtDNA),
  - female 2 / male 1 (chrX outside PAR),
  - female 0 / male 1 (chrY) — only male samples are included.
- **D.7.1 / D.7.2 `make_genDB_*`** Import GVCFs into GenomicsDB. Long intervals are imported as sub-batches (D.7.1); intervals shorter than the 2 Mb chunk size are grouped into short-segment sub-batches of up to 1000 segments each (D.7.2). This keeps GenomicsDB imports bounded in both memory and number of intervals per job.
- **D.8.1 / D.8.2 `GenotypeGVCFs_*`** Run `gatk GenotypeGVCFs` over each sub-batch.
- **D.9.1 `picardconcat`** If a (batch, ploidy) was chunked, concatenate the sub-batch GVCFs into one file per (batch, fploidy, mploidy).
- **D.9.2 `renameGVCF`** If there was only one sub-batch, rename it to the canonical output path instead.
- **D.10 `IndexGVCFs`** Index the final per-species, per-batch, per-ploidy GVCF.

## Design notes

- **Sentinel `done/` files** drive dependencies between targets instead of file timestamps. Sentinels are cheap to inspect and survive being moved across filesystems, which makes partial reruns straightforward.
- **Sharding and chunking** at three different layers (uBAM shards, calling batches, genotyping sub-batches) keep individual jobs short enough to fit cluster wall-time limits and let failed jobs be retried in isolation.
- **Per-region ploidy** is carried end-to-end through the regions table, the per-individual BED/intervals files, and the GenomicsDB workspace layout. This lets the same workflow handle autosomes, PAR, non-PAR X, Y, and mtDNA without per-species code paths.
- **Sex-aware calling** uses the genetic sex assigned from coverage (`gSEX` in `samples_coverage_stats.txt`, produced by `../annotate_sex_chr/`) rather than reported sex, so that mislabelled samples do not get the wrong ploidy.
- **Conditional fan-out** (`if os.path.isfile(...)`) means Block C and Block D only schedule jobs when their inputs already exist. Running `gwf run` repeatedly progressively adds new targets as upstream jobs complete.
