Scripts and data to recreate Figure 5.

### chain_to_beds.py
Takes a chain file (gzipped or not) and outputs three BED files:
  - <prefix>.unique.bed      : positions covered by exactly 1 chain
  - <prefix>.nonunique.bed   : positions covered by >1 chains
  - <prefix>.noalign.bed     : positions not covered by any chain

### callable_regions.py
For each species, loops over  BED files and computes total length per chromosome per bed and writes a table with columns:
  SPECIES, BATCH, FPLOIDY, MPLOIDY, BED_TYPE, CHR, N_CALLABLE
