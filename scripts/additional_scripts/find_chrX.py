# This script reads a <samples_file>  like subset_B_vasili.txt
# First, it reads the coverage for each contig longer than <min_chr_len> for each sample and normalizes these values per sample by the coverage of the longest contig
# Then, it iterates over reference genomes. For each reference, it identifies contigs that have coverage below 0.75 in at least one male. 
# These are candidate chrX contigs.
# Finally, these contigs are QC'ed: each of such contigs should have coverage below 0.75 i nmore than 50% of males and above 0.75 in more than 50% females.
# If a contig fails to meet at least one of the two conditions above (at least for F or M), it is not concidered a part of the X.
# If there are less than 50% of F or M with incorrect ploidy, the contig is believed to be chrX-derived while the samples are reported as potentially having wrong sex metadata.

# The script outputs an updated regions file with new ploidies for the discovered chrX contigs as well as for chrY and chrM if these are known and provided in the <special_contigs_file>.
# It also checks that the known chrX (if available for the reference) has the expected coverage.


import numpy as np
import pandas as pd
import argparse

#########################################################################################################################################
# Defining inputs and outputs
#########################################################################################################################################


# i use the argparse parser to feed in parameters

# These include 

# 1. the file with the special contigs which should look like this
# reference   contig  chrom
# Gorilla_gorilla_ssp NC_011120.1 chrM
# Gorilla_gorilla_ssp NC_073247.2 chrX
# Gorilla_gorilla_ssp NC_073248.2 chrY
# Macaca_mulatta_ssp  NC_005943.1 chrM
# Macaca_mulatta_ssp  NC_041774.1 chrX
# Macaca_mulatta_ssp  NW_021160381.1  chrX
# Macaca_mulatta_ssp  NW_021160382.1  chrX
# Macaca_mulatta_ssp  NW_021160383.1  chrX
# Macaca_mulatta_ssp  NW_021160384.1  chrX
# Macaca_mulatta_ssp  NC_027914.1 chrY
# Pan_paniscus_ssp    NC_001644.1 chrM
# Pan_paniscus_ssp    NC_073272.2 chrX
# Pan_paniscus_ssp    NC_073273.2 chrY
# Pan_troglodytes_ssp NC_001643.1 chrM
# Pan_troglodytes_ssp NC_072421.2 chrX
# Pan_troglodytes_ssp NC_072422.2 chrY
# Pithecia_pithecia_ssp   CM052647.1  chrX
# Plecturocebus_cupreus_ssp   JBDJOS010000622.1   chrM
# Plecturocebus_cupreus_ssp   CM080837.1  chrX
# Pongo_abelii_ssp    NC_002083.1 chrM
# Pongo_abelii_ssp    NC_072008.2 chrX
# Pongo_abelii_ssp    NC_072009.2 chrY
# Pongo_pygmaeus_ssp  NC_001646.1 chrM
# Pongo_pygmaeus_ssp  NC_072396.2 chrX
# Pongo_pygmaeus_ssp  NC_072397.2 chrY
# Theropithecus_gelada_ssp    NC_019802.1 chrM
# Theropithecus_gelada_ssp    NC_037689.1 chrX
# Trachypithecus_francoisi_ssp    NC_023970.1 chrM
# Trachypithecus_francoisi_ssp    NW_022681471.1  chrX


# 2. The subset file - this is the central file of out workflow. An example is subset_B_vasili.txt

# 3. The minimum length of a contig considered in the procedure. Default is 1Mb 


parser = argparse.ArgumentParser(description="")


parser.add_argument('-c', '--contigs', type=str, help='path to the file with special contigs (chrX, chrY and chrM)', 
	default = '~/special_contigs')

parser.add_argument('-s', '--subset', type=str, help='path to the file with the subset of samples', required=True)

parser.add_argument('-r', '--ref', type=str, help='the reference genome', required=True)

parser.add_argument('-m', '--minlen', type=float, help='minimum length of a contig to search for chrX derived contigs', 
    default = 1e6)


args = parser.parse_args()

special_contigs_file = args.contigs
samples_file = args.subset
ref = args.ref
min_chr_len = args.minlen


# There are also species and reference specufic file that follow a standard file naming. For now I keep them fixed but we can always make them as input parameters
cov_files = '~/{species}/cov/{id}.cov'
region_files = '~/{ref}/ref/regions_{ref}.txt'

# Outputs:
# the file below is the main output of the script. There are cases when it is not written for a given reference. That is a signal to read the logs!
updated_region_file =  '~/{ref}/ref/regions_{ref}_updated.txt'
# samples_out_file = '~/{ref}/ref/samples_with_potentially_wrong_sex.txt'
coverage_stats_out_file = '~/{ref}/ref/samples_coverage_stats.txt'


# this is a dummy variable which will be turned to False if the reported chrX has a non chrX coverange pattern
all_good = True

# to consider a contig haploid its relative coverage should be < coverage cutoff
cov_cutoff = 0.75

# to consider a contig Y-derived the fraction of its length covered should be < len_cutoff in females but >= len_cutoff in males
len_cutoff = 0.5

# a cutoff to decide if the contig is not diploid or some samples are wrongly sexed. This is the fraction of individuals having unexpected coverage
sample_fract = 0.5


#########################################################################################################################################
# Reading and processing input files
#########################################################################################################################################

print('Processing reference {ref}'.format(ref = ref))

# load the region file with F and M ploidy for each contig 
# this file will be modified based on a) known sex chromosomes and b) coverage in males and females
regions = pd.read_csv(region_files.format(ref = ref), sep = "\t")

original_regions = regions

# reading the file with special contigs (chrX, chrY and chrM) per reference
# this file has to be created manualy
special_contigs_df = pd.read_csv(special_contigs_file, sep = "\t", dtype = 'str')

special_contigs_df = special_contigs_df[special_contigs_df['reference'] == ref]

reported_X = special_contigs_df[special_contigs_df['chrom'] == 'chrX']['contig'].tolist()
reported_Y = special_contigs_df[special_contigs_df['chrom'] == 'chrY']['contig'].tolist()
reported_M = special_contigs_df[special_contigs_df['chrom'] == 'chrM']['contig'].tolist()


# reading the <subset> file and simplifying it to make a df with samples, their species, sex and reference genome
samples_df = pd.read_csv(samples_file, sep = "\t")
samples_df = samples_df[samples_df['REFERENCE_FOLDER'] == ref]
samples_df = samples_df[['IND_ID', 'FOLDER', 'SEX']].drop_duplicates()
samples_df = samples_df.reset_index()
samples_df = samples_df.drop(columns=['index'])


# counting individuals of each sex for each reference
sample_counts = samples_df.groupby(['SEX']).agg({'IND_ID': 'count'}).reset_index()
sample_counts = sample_counts.rename(columns={"SEX": "SEX", "IND_ID": "N"})


# I'm counting males and females per reference. If there is at least one of each I go on with estimating the candidate chrX and chrY contigs
# If one of the sexes is not represented I just update the regions file based on the provided info about sex chromosomes and mtDNA
N_males = sample_counts[(sample_counts['SEX'] == 'M')]['N'].tolist()
N_females = sample_counts[(sample_counts['SEX'] == 'F')]['N'].tolist()




#########################################################################################################################################
# Extracting coverage per contig per sample
#########################################################################################################################################


# reading individual cov files one by one, keeping only contigns longer than <min_chr_len>,
# extracting coverage for each such contig and normalizing coverage by the coverage of the longest contig
# we record absolute and relative coverage, contig and region names, and sample ID
# the mess with region names is to accomodate PARs and is perhaps redundant

print('Reading individual coverage files')

cov_list = []
cov_norm_list = []
len_list = []
sequenced_list = []
len_covered_list = []
contig_list = []
region_names_list = []
id_list = []
candidate_X = []
candidate_Y = []

bad_samples_df = pd.DataFrame()
bad_contigs_df = pd.DataFrame()


for sample_idx in range(samples_df.shape[0]):

    print(cov_files.format(species = samples_df['FOLDER'][sample_idx], id = samples_df['IND_ID'][sample_idx]))
    
    cov_df = pd.read_csv(cov_files.format(species = samples_df['FOLDER'][sample_idx], id = samples_df['IND_ID'][sample_idx]), 
                     sep = '\t',
                     names = ["chr_name", "chr", "start", "end", "covered", "sequenced", "cov"])


    cov_df['length'] = cov_df['end'] - cov_df['start']
    
    # keeping only contigs longer than the threshold
    cov_df = cov_df[cov_df['length'] > min_chr_len]

    length = cov_df['length'].tolist()
    len_list.extend(length)

    sequenced = cov_df['sequenced'].tolist()
    sequenced_list.extend(sequenced)

    len_covered = cov_df['covered'].tolist()
    len_covered_list.extend(len_covered)
    

    # using the coverage of the longest contig as the baseline to normalize
    # it can happen that the chrX is actually the longest one in the assembly. To prevent this from causing issues I make cov_df_tmp without chrX
    cov_df_tmp = cov_df[~cov_df['chr_name'].isin(reported_X)]
    baseline = cov_df_tmp[cov_df_tmp['length'] == max(cov_df_tmp['length'])]['cov'].tolist()[0]

    
    cov = cov_df['cov']
    cov_list.extend(cov)

    cov_norm_list.extend(round(cov / baseline, 3))

    contigs = cov_df['chr'].tolist()
    contig_list.extend(contigs)

    region_names = cov_df['chr_name'].tolist()
    region_names_list.extend(region_names)

    id_list.extend([samples_df['IND_ID'][sample_idx]] * len(contigs))


# baking it into a df
cov_per_chr_df = pd.DataFrame({'IND_ID': id_list,
              'contig': contig_list,
              'region_name': region_names_list,
              'cov': cov_list,
              'cov_norm': cov_norm_list,
              'length': len_list,
              'sequenced': sequenced_list,
              'len_covered_raw': len_covered_list})

# in the case there are more than one row for a given contig (so different regions like the PAR on the chrX) I will keep the coverage only for the 
# longest region of the contig
idx = cov_per_chr_df.groupby(['IND_ID', 'contig'])['length'].transform("max") == cov_per_chr_df['length']
cov_per_chr_df = cov_per_chr_df[idx]

cov_per_chr_df['len_covered'] = round(cov_per_chr_df['len_covered_raw'] / cov_per_chr_df['length'], 3)

cov_per_chr_df = pd.merge(cov_per_chr_df, samples_df, how="left", on=["IND_ID"])


print('writting the fraction of the genome covered per sample to file')

# # calculating the fraction of the entire genome covered by at least one read
# length_covered = cov_per_chr_df.groupby(['IND_ID', 'FOLDER']).apply( lambda df,a,b: sum(df[a] * df[b]) / sum(df[b]), a = 'len_covered', b = 'length').reset_index()
# length_covered.to_csv(length_covered_file.format(ref = ref), index=False, sep = '\t', header=False)

print('\n')
    
if len(N_males) > 0 and len(N_females) > 0:

    print('There are {NF} females and {NM} males for this reference. Proceeding with chrX inference.'.format(NF = N_females[0], NM = N_males[0]))

    #########################################################################################################################################
    # Finding contigs that can potentially be chrX-derived
    #########################################################################################################################################

    # get a list of contigs which potentially come from the X chromosome 
    # the conditions is that they have less than 0.75 coverage in at least one male
    candidate_X = cov_per_chr_df[(cov_per_chr_df['cov_norm'] < cov_cutoff) & (cov_per_chr_df['SEX'] == 'M')][['contig']].drop_duplicates()
    candidate_X = candidate_X['contig'].tolist()
    
    # a df with coverage of those contings in each sample
    candidate_X_df = cov_per_chr_df[(cov_per_chr_df['contig'].isin(candidate_X))]
    
    # find the candidate chrX contigs that don't have expected coverage in at least one sample
    flaged_contigs = candidate_X_df[(candidate_X_df['cov_norm'] < cov_cutoff) & (candidate_X_df['SEX'] == 'F') | (candidate_X_df['cov_norm'] >= cov_cutoff) & (candidate_X_df['SEX'] == 'M')]
    
    # simplifying the df so that it only has the name of the candidate chrX contig and the corresponding reference name
    candidate_X_df = candidate_X_df[['contig', 'region_name']].drop_duplicates()


    #########################################################################################################################################
    # Checking the candidate contigs. If a contig is not meeting the criteria in most of the samples it is removed (its ploidy will not be changed)
    # If a contig is suspecious in only few samples it will still be treated as chrX but the problematic cases will be reported
    #########################################################################################################################################


    # check if there are any problematic contigs
    
    if flaged_contigs.shape[0] == 0:
        print('All the candidate chrX contigs have correct coverage in the all samples')
    else:
        # for each problematic contig we count the number of individuals, F and M separately, that have incorrect coverage 
        flaged_contigs_counts = flaged_contigs.groupby(['contig', 'region_name', 'SEX']).agg({'IND_ID': 'count'}).reset_index()
        
        # and add the total count of each sex for each reference
        flaged_contigs_counts = pd.merge(flaged_contigs_counts, sample_counts, how="left", on=[ 'SEX'])
    
        # if the count of individuals with incorrect coverage for a given contig is below 50% we conclude that the contig is ok but some samples are mis-labeled
        # we record these samples in a separate df and write them to file
        bad_samples_df = flaged_contigs_counts[flaged_contigs_counts['IND_ID'] / flaged_contigs_counts['N'] < sample_fract]
        bad_contigs_df = flaged_contigs_counts[flaged_contigs_counts['IND_ID'] / flaged_contigs_counts['N'] >= sample_fract]
        
        # if bad_samples_df.shape[0] > 0:
        #     print('Some samples have potentially wrong sex metadata. \nThese are written to ', samples_out_file.format(ref = ref))
    
        #     bad_samples_df = flaged_contigs[(flaged_contigs['contig'].isin(bad_samples_df['contig'])) & (flaged_contigs['region_name'].isin(bad_samples_df['region_name']))]
        #     bad_samples_df.to_csv(samples_out_file.format(ref = ref), index=False, sep = '\t')
    
        if bad_contigs_df.shape[0] > 0:
            # remove the contigs with unexpected coverage in more than or equal to 50% of samples by sex from the list of candidate contigs
            # candidate_X_df = candidate_X_df[~((candidate_X_df['contig'].isin(bad_contigs_df['contig'])) & (candidate_X_df['region_name'].isin(bad_contigs_df['region_name'])))] 
            bad_candidate_X = bad_contigs_df['contig'].tolist()
            candidate_X = list(set(candidate_X) - set(bad_candidate_X))
            
            # As the original criteria for being candidate_X is cov < 0.75 in at least one male the 'bad_contigs' may include Y-chr
            candidate_Y = cov_per_chr_df[cov_per_chr_df['contig'].isin(bad_contigs_df[bad_contigs_df['SEX'] == 'F']['contig'])]
            if candidate_Y.shape[0] > 0:
                candidate_Y = candidate_Y.groupby(['contig', 'SEX']).agg({'cov_norm': 'mean', 'len_covered': 'mean'}).reset_index()
                candidate_Y = candidate_Y.pivot(index='contig', columns='SEX', values='len_covered')
                candidate_Y = candidate_Y[(candidate_Y['F'] < len_cutoff) & (candidate_Y['M'] >= len_cutoff)].reset_index()['contig'].tolist()
                bad_candidate_Y = bad_contigs_df[(bad_contigs_df['contig'].isin(candidate_Y)) & bad_contigs_df['SEX'] == 'M']['contig'].tolist()
                candidate_Y = list(set(candidate_Y) - set(bad_candidate_Y))
                
                reported_Y = set(list(set(candidate_Y)) + list(set(reported_Y))) 

                bad_candidate_X = list(set(bad_candidate_X) - set(reported_Y))    

            if len(bad_candidate_X) > 0:
                print('Some contigs have coverage below 0.75 in some males but do not meet all the criteria.')
                print(bad_candidate_X)
                


    reported_X = set(list(set(candidate_X)) + list(set(reported_X)))

    reported_X = list(set(reported_X) - set(reported_Y))

    


    # check the final list of chrX and chrY contigs. 
    # chrX should have mean coverage >= 0.75 in F and < 0.75 in males
    # chrY should have mean coverage (across samples) < 0.75 in both males and females
    # if any contig is not meeting these criteria, there will be a warning and the region file will not be updated
    def get_check_df(candidates):
        check_df = cov_per_chr_df[cov_per_chr_df['contig'].isin(candidates)].groupby(['contig', 'SEX']).agg({'cov_norm': 'mean'}).reset_index()
        check_df = check_df.pivot(index='contig', columns='SEX', values='cov_norm')
        return check_df

    if len(reported_X) > 0:
        chrX_check = get_check_df(reported_X)
        chrX_check = chrX_check[~((chrX_check['F'] >= cov_cutoff) & (chrX_check['M'] < cov_cutoff))]

        if chrX_check.shape[0] > 0:
            print('Some contigs reported as chrX have an unmatching coverage profile')
            print(chrX_check.to_string())
            all_good = False

    # if len(reported_Y) > 0:
    #     chrY_check = get_check_df(reported_Y)
    #     chrY_check = chrY_check[~((chrY_check['F'] < 0.75) & (chrY_check['M'] < 0.75))]

    #     if chrY_check.shape[0] > 0:
    #         print('Some contigs reported as chrY have an unmatching coverage profile')
    #         print(chrY_check.to_string())
    #         all_good = False


# updating the region files    
    
if len(reported_X) > 0:
    print('There are {N} contigs reported / identified as coming from the chrX. Setting their ploidy to 1 in males'.format(N = len(reported_X)))
    # regions.loc[(regions['chrom'].isin(reported_X)), 'male_ploidy'] = 1
    regions.loc[(regions['chrom'].isin(reported_X)) & (regions['region'] != 'PAR1') & (regions['region'] != 'PAR2'), 'male_ploidy'] = 1
else:
    print('There is no reported / identified chrX for this reference')


if len(reported_Y) > 0:
    print('There are {N} contigs reported / identified as coming from the chrY. Setting their ploidy to 1 in males and 0 in females'.format(N = len(reported_Y)))
    regions.loc[regions['chrom'].isin(reported_Y), 'male_ploidy'] = 1
    regions.loc[regions['chrom'].isin(reported_Y), 'female_ploidy'] = 0
else:
    print('There is no reported / identified chrY for this reference')


if len(reported_M) > 0:
    print('There are {N} contigs reported as coming from the chrM. Setting their ploidy to 1 in males and 1 in females'.format(N = len(reported_M)))
    regions.loc[regions['chrom'].isin(reported_M), 'male_ploidy'] = 1
    regions.loc[regions['chrom'].isin(reported_M), 'female_ploidy'] = 1
else:
    print('There is no reported chrM for this reference')
    

if all_good:
    print('Writting the updated regions file')
    if ref in ['Gorilla_gorilla_ssp', 'Homo_sapiens_ssp', 'Pan_paniscus_ssp', 'Pan_troglodytes_ssp', 'Pongo_abelii_ssp', 'Pongo_pygmaeus_ssp', 'Symphalangus_syndactylus_ssp']:
        original_regions.to_csv(updated_region_file.format(ref = ref), index=False, sep = '\t')
    else:
        regions.to_csv(updated_region_file.format(ref = ref), index=False, sep = '\t')
else:
    print('There is some inconsistency between the provided sex chromosome contig names and the coverage of males and females.\n No updated regions file written.')



# caclulating coverage statistics for the autosomes, chrX and chrY

# locating the corresponding contigs
nonautosomes = regions[regions['male_ploidy'] != 2]['chrom'].tolist()
chrX = regions[(regions['female_ploidy'] == 2) & (regions['male_ploidy'] == 1)]['chrom'].tolist()
chrY = regions[(regions['female_ploidy'] == 0) & (regions['male_ploidy'] == 1)]['chrom'].tolist()


# subsetting cov_per_chr_df
autosomes_df = cov_per_chr_df[~cov_per_chr_df['contig'].isin(nonautosomes)]
chrX_df = cov_per_chr_df[cov_per_chr_df['contig'].isin(chrX)]
chrY_df = cov_per_chr_df[cov_per_chr_df['contig'].isin(chrY)]


# a function to summarise coverage for an input dataframe and add a postfix to its columns for downstream merging

def summarise_cov(df, postfix):
    df = df.groupby(['IND_ID']).agg({'length': 'sum', 'sequenced': 'sum', 'len_covered_raw': 'sum'}).reset_index()
    df['cov'] = round(df['sequenced'] / df['len_covered_raw'], 3)
    df['cov_len'] = round(df['len_covered_raw'] / df['length'], 3)
    df = df.rename(columns={c: c + '_' + postfix for c in df.columns if c not in ['IND_ID']})
    return df

# summarise over autosomes and add the sex from metadata

summarise_cov_df = summarise_cov(autosomes_df, 'A')
summarise_cov_df = pd.merge(summarise_cov_df, samples_df[['IND_ID', 'SEX']], how="left", on=[ 'IND_ID'])

# if chr X contigs exist, add cov info on those and genetically inferred sex
if chrX_df.shape[0] > 0:
    summarise_cov_df = pd.merge(summarise_cov_df, summarise_cov(chrX_df, 'X'), how="left", on=[ 'IND_ID'] )
    summarise_cov_df['rel_X_cov'] = round(summarise_cov_df['cov_X'] / summarise_cov_df['cov_A'], 3)
    summarise_cov_df['gSEX'] = ['M' if c < cov_cutoff else 'F' for c in summarise_cov_df['rel_X_cov'] ] 
else:
    summarise_cov_df['gSEX'] = summarise_cov_df['SEX']


# if chr Y contigs exist, add cov info on those
if chrY_df.shape[0] > 0:
    summarise_cov_df = pd.merge(summarise_cov_df, summarise_cov(chrY_df, 'Y'), how="left", on=['IND_ID'])


summarise_cov_df.to_csv(coverage_stats_out_file.format(ref = ref), index=False, sep = '\t')
