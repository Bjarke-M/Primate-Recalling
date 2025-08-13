import cyvcf2
import numpy as np
import allel
import pandas as pd
from collections import defaultdict
import sys

def read_ploidy_info(ploidy_file):
    """Read ploidy information and identify sex chromosomes"""
    regions = pd.read_csv(ploidy_file, sep='\t')
    
    # Identify X chromosomes
    X_chroms = regions.loc[(regions['female_ploidy'] == 2) & (regions['male_ploidy'] == 1), 'chrom'].values
    
    # Identify Y chromosomes
    Y_chroms = regions.loc[(regions['female_ploidy'] == 0) & (regions['male_ploidy'] == 1), 'chrom'].values
    
    # All other chromosomes are autosomes
    all_chroms = regions['chrom'].values
    autosomes = [c for c in all_chroms if c not in X_chroms and c not in Y_chroms]
    
    return autosomes, X_chroms, Y_chroms

def read_coverage_info(mask_file):
    """Read coverage information and filter for samples in BCF"""
    mask=pd.read_csv(mask_file, sep= '\t', names = ['Chromosome', 'Start', 'End'])
    # Filter to only include samples present in the BCF file
    mask['interval_length'] = mask['End'] - mask['Start']
    total_callable_length = mask.groupby('Chromosome')['interval_length'].sum().reset_index(name='Length')
    return total_callable_length

def count_heterozygous_sites(vcf_path):
    """Count heterozygous sites per individual per chromosome"""
    
    # Load VCF file
    vcf = cyvcf2.VCF(vcf_path)
    samples = vcf.samples
    
    # Dictionary to store data by chromosome
    chr_data = defaultdict(lambda: {'variants': [], 'genotypes': []})
    
    # Read variants and group by chromosome
    for variant in vcf:
        chrom = variant.CHROM
        
        # Collect variant info
        chr_data[chrom]['variants'].append({
            'CHROM': variant.CHROM,
            'POS': variant.POS,
            'REF': variant.REF,
            'ALT': ','.join(variant.ALT) if variant.ALT else '.',
            'QUAL': variant.QUAL,
            'ID': variant.ID
        })
        
        # Get genotype data (just the alleles, not the phasing)
        gt = variant.genotype.array()
        chr_data[chrom]['genotypes'].append(gt[:, :2])
    
    # Initialize results dictionary
    het_counts = {sample: {} for sample in samples}
    
    # Process each chromosome
    for chrom, data in chr_data.items():
        if not data['genotypes']:  # Skip empty chromosomes
            continue
            
        # Convert to numpy array and create GenotypeArray
        genotype_array = np.array(data['genotypes'])
        gt_array = allel.GenotypeArray(genotype_array)
        
        # Count variants on this chromosome before filtering
        n_sites_total = gt_array.shape[0]
        
        # Count alleles to identify polymorphic sites
        ac = gt_array.count_alleles()
        # A site is Biallelic if it has more 2 allels
        is_biallelic = ac.allelism() == 2
        
        # Combine filters: keep sites that pass filter
        mask = is_biallelic
        
        # Apply filters
        gt_array_filtered = gt_array[mask]
        
        # Count variants after filtering
        n_variants = gt_array_filtered.shape[0]
        
        if n_variants == 0:  # Skip if no variants pass filters
            for i, sample in enumerate(samples):
                het_counts[sample][chrom] = 0
            continue
        
        # Check which genotypes are heterozygous
        is_het = gt_array_filtered.is_het()
        
        # Count heterozygous sites per sample (sum across variants axis=0)
        n_het_per_sample = np.sum(is_het, axis=0)
        
        # Store results for each sample
        for i, sample in enumerate(samples):
            het_counts[sample][chrom] = int(n_het_per_sample[i])
    return het_counts, samples

def get_sexes(samples, sex_file):
    sexes = pd.read_csv(sex_file, sep='\t')
    sexes = sexes[['BIOSAMPLE_ID', 'GENETIC_SEX', 'GVCF_FOLDER', 'PDGP_ID']]
    # Filter for samples present in either 'Sample' or 'PDGP_ID'
    mask = sexes['BIOSAMPLE_ID'].isin(samples) | sexes['PDGP_ID'].isin(samples)
    return sexes[mask]

###### load in files #####
bcf_file = sys.argv[1]
ploidy_file = sys.argv[2]
mask_file = sys.argv[3]
sexes_file= sys.argv[4]
out = sys.argv[5]



autosomes, X_chroms, Y_chroms = read_ploidy_info(ploidy_file)


het_counts, samples = count_heterozygous_sites(bcf_file)

total_callable = read_coverage_info(mask_file)

sexes = get_sexes(samples,sexes_file)


# Prepare a list to collect results
results = []

# Convert total_callable to a dictionary for quick lookup
callable_dict = dict(zip(total_callable['Chromosome'], total_callable['Length']))

# For each sample, process chromosomes
for sample in samples:
    # Get the genetic sex and BIOSAMPLE_ID for this sample (using BIOSAMPLE_ID or PDGP_ID)
    sex_row = sexes[(sexes['BIOSAMPLE_ID'] == sample) | (sexes['PDGP_ID'] == sample)]
    if not sex_row.empty:
        genetic_sex = sex_row.iloc[0]['GENETIC_SEX']
        biosample_id = sex_row.iloc[0]['BIOSAMPLE_ID']
        pdgp_id = sex_row.iloc[0]['PDGP_ID']
        species = sex_row.iloc[0]['GVCF_FOLDER'] 
    else:
        genetic_sex = None  # If not found, treat as unknown
        biosample_id = sex_row.iloc[0]['BIOSAMPLE_ID']
        pdgp_id = sex_row.iloc[0]['PDGP_ID']

    for chrom in het_counts[sample]:
        # Skip X and Y chromosomes for males
        if genetic_sex == 'M' and (chrom in X_chroms or chrom in Y_chroms):
            continue

        het_count = het_counts[sample][chrom]
        callable_length = callable_dict.get(chrom, np.nan)
        het_per_bp = het_count / callable_length if callable_length and callable_length > 0 else np.nan
        if chrom in autosomes:
            chrom_type = 'A'
        elif chrom in X_chroms:
            chrom_type = 'X'
        elif chrom in Y_chroms:
            chrom_type = 'Y'
        else:
            chrom_type = 'Unknown'

        results.append({
            'Sample': pdgp_id,
            'BIOSAMPLE_ID': biosample_id,
            'Chromosome': chrom,
            'Chrom_Type': chrom_type,
            'Het_Count': het_count,
            'Callable_Length': callable_length,
            'Genetic_sex': genetic_sex,
            'Het_per_bp': het_per_bp,
            'Species': species
        })

# Create a DataFrame from the results
df_results = pd.DataFrame(results)

df_results.to_csv(out, index=False)

