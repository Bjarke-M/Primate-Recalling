# figure5

Callable regions from chain-file alignments (Figure 5).

## Scripts

- `chain_to_beds.py` — takes a chain file (gzipped or plain) and writes three BED files:
  - `<prefix>.unique.bed` — positions covered by exactly one chain.
  - `<prefix>.nonunique.bed` — positions covered by more than one chain.
  - `<prefix>.noalign.bed` — positions not covered by any chain.
- `callable_regions.py` — for each species, sums BED interval lengths per chromosome and per BED type. Output columns: `SPECIES, BATCH, FPLOIDY, MPLOIDY, BED_TYPE, CHR, N_CALLABLE`.
- `fig5.R` — plots Figure 5.

## Data

- `../data/alnRegionHet.zip` — per-species outputs of `chain_to_beds.py` and `callable_regions.py`.
