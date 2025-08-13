import sys 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as patches


##################### Helper functions #####################
def mask_coverage_fast(df, mask):
    df = df.copy()
    df['total_overlap'] = 0
    
    for chr_name in df['Chromosome'].unique():
        # Get intervals for this chromosome
        roh_chr = df[df['Chromosome'] == chr_name]
        mask_chr = mask[mask['Chromosome'] == chr_name]
        
        if mask_chr.empty:
            continue
            
        # Create IntervalIndex for mask intervals
        mask_intervals = pd.IntervalIndex.from_arrays(
            mask_chr['Start'], mask_chr['End'], closed='both'
        )
        
        # For each ROH interval, find overlapping mask intervals and calculate overlap
        for idx in roh_chr.index:
            roh_start = roh_chr.loc[idx, 'Start']
            roh_end = roh_chr.loc[idx, 'End']
            
            # Find overlapping mask intervals
            overlapping_idx = mask_intervals.overlaps(pd.Interval(roh_start, roh_end, closed='both'))
            overlapping_masks = mask_chr[overlapping_idx]
            
            if not overlapping_masks.empty:
                # Vectorized overlap calculation
                overlap_starts = np.maximum(roh_start, overlapping_masks['Start'])
                overlap_ends = np.minimum(roh_end, overlapping_masks['End'])
                total_overlap = (overlap_ends - overlap_starts).sum()
                df.loc[idx, 'total_overlap'] = total_overlap
    
    return df

def create_id_mapping(sex_df):
    """
    Create a mapping from any ID format to BIOSAMPLE ID (canonical)
    """
    id_mapping = {}
    
    # Add BIOSAMPLE IDs (map to themselves)
    for biosample in sex_df['Sample'].unique():
        if pd.notna(biosample):
            id_mapping[biosample] = biosample
    
    # Add PD ID mappings
    if 'PDGP_ID' in sex_df.columns:
        for _, row in sex_df.iterrows():
            if pd.notna(row['PDGP_ID']) and pd.notna(row['Sample']):
                # Handle different PD formats
                pd_id = row['PDGP_ID']
                biosample_id = row['Sample']
                
                # Map different PD formats to BIOSAMPLE
                id_mapping[pd_id] = biosample_id
                id_mapping[pd_id.replace('_', ' ')] = biosample_id  # "PD_0282" -> "PD 0282"
                id_mapping[pd_id.replace(' ', '_')] = biosample_id  # "PD 0282" -> "PD_0282"
    
    return id_mapping

def standardize_sample_ids(df, id_mapping):
    """
    Convert all sample IDs to their canonical BIOSAMPLE form
    """
    df = df.copy()
    
    # Try to map each sample ID
    def map_sample(sample):
        # Clean up the sample ID
        sample_str = str(sample).strip()
        
        # Try direct mapping
        if sample_str in id_mapping:
            return id_mapping[sample_str]
        
        # Try with space/underscore variations
        sample_underscore = sample_str.replace(' ', '_')
        if sample_underscore in id_mapping:
            return id_mapping[sample_underscore]
            
        sample_space = sample_str.replace('_', ' ')
        if sample_space in id_mapping:
            return id_mapping[sample_space]
        
        # If no mapping found, return original
        return sample_str
    
    df['Sample'] = df['Sample'].apply(map_sample)
    return df

##################### load in files #####################
ROH_file = sys.argv[1]
mask_file = sys.argv[2]
regions_file= sys.argv[3]
sex_file = sys.argv[4]
species_ssp = sys.argv[5]
outfile = sys.argv[6]
out_fig_seg= sys.argv[7]
out_fig_dist = sys.argv[8]

#ROH file from plink
ROH = pd.read_csv(ROH_file, sep='\s+')

# Check if FID and IID columns are identical
if (ROH['FID'] == ROH['IID']).all():
    ROH = ROH.rename(columns={'FID':'Sample','CHR':'Chromosome','POS1':'Start','POS2':'End','KB':'Length'})
else:
    ROH['Sample'] = ROH['FID'].astype(str) + '_0' + ROH['IID'].astype(str)
    ROH = ROH.rename(columns={'CHR':'Chromosome','POS1':'Start','POS2':'End','KB':'Length'})

ROH['Length']=ROH['Length']*1e3 #from kb to b
ROH = ROH[['Sample','Chromosome','Start','End','Length']]

# Mask file
mask = pd.read_csv(mask_file, sep= '\t', names = ['Chromosome', 'Start', 'End'])

# Regions file to get sex chromosomes
regions = pd.read_csv(regions_file, sep='\t')

# Sex file, to get genetic sex
sexes = pd.read_csv(sex_file, sep='\t')
sexes = sexes.loc[sexes['GVCF_FOLDER'] == species_ssp, ['BIOSAMPLE_ID', 'GENETIC_SEX', 'GVCF_FOLDER', 'PDGP_ID']]
sexes = sexes.rename(columns={'BIOSAMPLE_ID':'Sample'})

# Create ID mapping
id_mapping = create_id_mapping(sexes)
print("ID mapping created:", id_mapping)

# Standardize ROH sample IDs to BIOSAMPLE format
ROH = standardize_sample_ids(ROH, id_mapping)

# Get all unique BIOSAMPLE IDs (canonical IDs only)
all_samples = sorted(sexes['Sample'].unique())
all_chromosomes = sorted(mask['Chromosome'].unique())

print("All samples (BIOSAMPLE IDs):", all_samples)

##################### Remove ROH with len(ROH) < 1mb #####################
ROH = ROH[ROH['Length'] >= 1000000]

##################### Combine mask and ROH #####################
df = mask_coverage_fast(ROH, mask)

##################### Get x and y chromosomes #####################
# Identify X chromosomes
X_chroms = regions.loc[(regions['female_ploidy'] == 2) & (regions['male_ploidy'] == 1), 'chrom'].values
# Identify Y chromosomes
Y_chroms = regions.loc[(regions['female_ploidy'] == 0) & (regions['male_ploidy'] == 1), 'chrom'].values

print("X chromosomes:", X_chroms)
print("Y chromosomes:", Y_chroms)

##################### Merge with sex information #####################
# Simple merge since sample IDs are now standardized
dataframe = pd.merge(df, sexes, on='Sample', how='left')

# Remove X and Y chromosome in males
dataframe = dataframe.loc[
    ~((dataframe['Chromosome'].isin(X_chroms) & (dataframe['GENETIC_SEX'] == 'M')) |
        (dataframe['Chromosome'].isin(Y_chroms)))]

##################### Calculate fROH #####################
dataframe.loc[:, 'cov_frac'] = dataframe['total_overlap'] / dataframe['Length']

##################### Length distribution of RoHs #####################
froh = dataframe.loc[dataframe['cov_frac'] > 0.5]

# Get unique samples that passed the filter
n_samples = len(all_samples)
n_cols = 3
n_rows = (n_samples + n_cols - 1) // n_cols

# Create figure and subplots
fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
if n_rows * n_cols > 1:
    axes = axes.flatten()
else:
    axes = [axes]

# Plot histogram for each sample
for idx, sample in enumerate(all_samples):
    sample_data = froh[froh['Sample'] == sample]['Length']
    
    if len(sample_data) > 0:
        axes[idx].hist(sample_data, bins=30, edgecolor='black')
    else:
        axes[idx].text(0.5, 0.5, 'No ROHs passing filter', 
                      ha='center', va='center', transform=axes[idx].transAxes)
        axes[idx].set_xlim(0, 1)
        axes[idx].set_ylim(0, 1)
    
    axes[idx].set_title(f'Sample: {sample}')
    axes[idx].set_xlabel('Length (bp)')
    axes[idx].set_ylabel('Count')

# Remove empty subplots if any
for idx in range(len(all_samples), len(axes)):
    fig.delaxes(axes[idx])

plt.tight_layout()
plt.savefig(out_fig_dist, format="pdf")  
plt.close()

##################### Segment plot #####################
df = dataframe.loc[dataframe['cov_frac'] > 0.5]

# Use the complete lists of samples and chromosomes
unique_samples = all_samples
unique_chromosomes = all_chromosomes

# Calculate dynamic figure size based on data
base_width_per_chrom = 40 / 24  # ~1.67 units per chromosome
base_height_per_sample = 8 / 18  # ~0.44 units per sample

# Calculate figure dimensions with minimum sizes
fig_width = max(12, len(unique_chromosomes) * base_width_per_chrom)
fig_height = max(4, len(unique_samples) * base_height_per_sample + 2)

# Create figure with dynamic proportions
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Define colors for different samples
colors_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
sample_colors = {sample: colors_list[i % len(colors_list)] 
                 for i, sample in enumerate(unique_samples)}

# Calculate chromosome boundaries
chromosome_info = {}
for chrom in unique_chromosomes:
    mask_chrom = mask[mask['Chromosome'] == chrom]
    if not mask_chrom.empty:
        chromosome_info[chrom] = {
            'min': 0,
            'max': mask_chrom['End'].max() + 5000000
        }
    else:
        chromosome_info[chrom] = {
            'min': 0,
            'max': 250000000  # Default
        }

# Calculate cumulative positions for chromosomes
chromosome_offsets = {}
cumulative_pos = 0
chromosome_spacing = 2000000

for chrom in unique_chromosomes:
    chromosome_offsets[chrom] = cumulative_pos
    cumulative_pos += chromosome_info[chrom]['max'] + chromosome_spacing

# Plot the segments
segment_height = 0.7
sample_y_positions = {sample: i for i, sample in enumerate(unique_samples)}

# Draw chromosome backgrounds first
for chrom in unique_chromosomes:
    offset = chromosome_offsets[chrom]
    chrom_length = chromosome_info[chrom]['max']
    
    for sample in unique_samples:
        y_pos = sample_y_positions[sample]
        rect = patches.Rectangle((offset, y_pos - segment_height/2), 
                               chrom_length, segment_height,
                               linewidth=0.5, 
                               edgecolor='gray',
                               facecolor='lightgray',
                               alpha=0.2,
                               zorder=1)
        ax.add_patch(rect)

# Draw the actual segments
for _, row in df.iterrows():
    sample = row['Sample']
    chrom = row['Chromosome']
    start = row['Start'] + chromosome_offsets[chrom]
    end = row['End'] + chromosome_offsets[chrom]
    
    if sample in sample_y_positions:
        y_pos = sample_y_positions[sample]
        rect = patches.Rectangle((start, y_pos - segment_height/2), 
                               end - start, segment_height,
                               linewidth=0, 
                               edgecolor='none',
                               facecolor=sample_colors[sample],
                               alpha=0.8,
                               zorder=2)
        ax.add_patch(rect)

# Add chromosome labels
label_fontsize = min(10, max(6, 240 / len(unique_chromosomes)))
for chrom in unique_chromosomes:
    offset = chromosome_offsets[chrom]
    chrom_length = chromosome_info[chrom]['max']
    
    ax.text(offset + chrom_length/2, len(unique_samples) + 0.5, 
            chrom.replace('NC_', '').replace('.2', '').replace('.1', ''), 
            ha='center', va='bottom', fontsize=label_fontsize, fontweight='bold')
    
    if chrom != unique_chromosomes[0]:
        ax.axvline(x=offset - chromosome_spacing/2, 
                  color='darkgray', linestyle='--', alpha=0.5, linewidth=1)

# Set axes
ytick_fontsize = min(10, max(6, 180 / len(unique_samples)))
ax.set_yticks(list(sample_y_positions.values()))
ax.set_yticklabels(list(sample_y_positions.keys()), fontsize=ytick_fontsize)

ax.set_xlabel('Genomic Position (Mb)', fontsize=12)
ax.set_ylabel('Samples', fontsize=12)

ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.0f}'))

ax.set_title('Homozygous/Autozygous Regions Across Chromosomes', 
             fontsize=16, pad=20, fontweight='bold')

ax.set_xlim(-chromosome_spacing, cumulative_pos - chromosome_spacing)
ax.set_ylim(-1, len(unique_samples) + 1)

# Styling
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, axis='x', alpha=0.2, linestyle='-', linewidth=0.5)
ax.set_axisbelow(True)

# Legend
legend_cols = min(3, max(1, len(unique_samples) // 6))
max_legend_items = min(10, len(unique_samples))

legend_patches = [patches.Patch(color=sample_colors[sample], label=sample, alpha=0.8) 
                  for sample in unique_samples[:max_legend_items]]
if len(unique_samples) > max_legend_items:
    legend_patches.append(patches.Patch(color='white', label=f'... (+{len(unique_samples) - max_legend_items} more)'))
    
ax.legend(handles=legend_patches, loc='upper right', frameon=True, 
          fancybox=True, shadow=False, ncol=legend_cols, fontsize=9)

plt.tight_layout()
plt.savefig(out_fig_seg, format='pdf', dpi=300, bbox_inches='tight')
plt.close()

##################### Summary dataframe #####################
# Callable length of each chromosome
mask['interval_length'] = mask['End'] - mask['Start']
total_callable_length = mask.groupby('Chromosome')['interval_length'].sum().reset_index(name='Length')
mask = mask.drop(columns=['interval_length'])  # Clean up temporary column

# Create complete sample-chromosome combinations
all_combinations = pd.DataFrame([(s, c) for s in all_samples for c in all_chromosomes], 
                                columns=['Sample', 'Chromosome'])

all_combinations['Chrom_Type'] = all_combinations['Chromosome'].apply(
    lambda x: 'X' if x in X_chroms else ('Y' if x in Y_chroms else 'A')
)

sex_mapping = sexes[['Sample', 'GENETIC_SEX']].drop_duplicates()
all_combinations = pd.merge(all_combinations, sex_mapping, on='Sample', how='left')

# Rename the column if you want it to be just 'Sex' instead of 'GENETIC_SEX'
all_combinations = all_combinations.rename(columns={'GENETIC_SEX': 'Sex'})


########### Remove X in Males ###########
all_combinations = all_combinations.loc[
    ~(all_combinations['Chromosome'].isin(X_chroms) & (all_combinations['Sex'] == 'M'))]

all_combinations['species'] = species_ssp

# Process filtered data
filtered_data = dataframe.loc[dataframe['cov_frac'] > 0.5]

if not filtered_data.empty:
    grouped_df = filtered_data.groupby(['Sample', 'Chromosome'])
    length_sum = grouped_df['Length'].sum().reset_index(name='total_length')
    cov_frac_median = grouped_df['cov_frac'].median().reset_index(name='median_coverage')
    count_per_group = grouped_df.size().reset_index(name='count')
else:
    length_sum = pd.DataFrame(columns=['Sample', 'Chromosome', 'total_length'])
    cov_frac_median = pd.DataFrame(columns=['Sample', 'Chromosome', 'median_coverage'])
    count_per_group = pd.DataFrame(columns=['Sample', 'Chromosome', 'count'])

# Build summary dataframe
summary_df = all_combinations.copy()

# Simple merges since all IDs are now standardized
summary_df = pd.merge(summary_df, length_sum, on=['Sample', 'Chromosome'], how='left')
summary_df['total_length'] = summary_df['total_length'].fillna(0)

summary_df = pd.merge(summary_df, total_callable_length, on='Chromosome', how='left')
summary_df = summary_df.rename(columns={'total_length': 'total_roh_len', 'Length': 'callable_len'})

summary_df['froh'] = summary_df['total_roh_len'] / summary_df['callable_len']
summary_df['froh'] = summary_df['froh'].fillna(0)

summary_df = pd.merge(summary_df, cov_frac_median, on=['Sample', 'Chromosome'], how='left')
summary_df['median_coverage'] = summary_df['median_coverage'].fillna(0)

summary_df = pd.merge(summary_df, count_per_group, on=['Sample', 'Chromosome'], how='left')
summary_df['count'] = summary_df['count'].fillna(0).astype(int)

# Save summary
summary_df.to_csv(outfile, index=False)
