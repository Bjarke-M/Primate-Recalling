import os
import sys
import pandas as pd
import subprocess

#### load data
length_of_chunk = int(sys.argv[1])
group           = sys.argv[2]

regions = pd.read_table(sys.argv[3])
regions = regions.sort_values(by=["batch", "chrom", "start"])
batches = list(regions.batch)
female_ploidy = list(regions.female_ploidy)
male_ploidy = list(regions.male_ploidy)
regions_list = list(regions.region)
chrom_list = list(regions.chrom)
start_list = [jj + 1 for jj in list(regions.start)]
end_list = list(regions.end)

inds        = sys.argv[4].split(",")
sexes       = sys.argv[5].split(",")
ref_folder  = sys.argv[6]
gvcf_folder = sys.argv[6]

#### delete old folder and make new floder
subprocess.run(["rm", "-fr", f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/geno_beds_and_intervals/"])
subprocess.run(["mkdir", "-p", f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/geno_beds_and_intervals/"])

#### make metadata files for individual calling
# for j in range(len(inds)):
#     #### make set_of_ploidy file for each individual based on sex
#     if sexes[j] == "F":
#         ploidies_list = list(regions.female_ploidy)
#     else:
#         ploidies_list = list(regions.male_ploidy)
#     for batch in sorted(list(set(regions.batch))):
#         iis = [index for index in range(len(batches)) if batches[index] == batch]
#         set_of_ploidies = sorted(list(set([ploidies_list[ii] for ii in iis if ploidies_list[ii] != 0])))
    
#         ## get ploidy of regions based on sex of individual
#         for pl in set_of_ploidies:
#             chromosomes = [chrom_list[ii] for ii in iis if ploidies_list[ii] == pl]
#             starts      = [start_list[ii] for ii in iis if ploidies_list[ii] == pl]
#             ends        = [end_list[ii] for ii in iis if ploidies_list[ii] == pl]

#             oo = open(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/geno_beds_and_intervals/{inds[j]}_batch_{batch}_ploidy_{pl}.bed", "w")
#             oo.write("\n".join(chromosomes[jj] + "\t" + str(starts[jj] - 1) + "\t" + str(ends[jj]) for jj in range(len(chromosomes))))
#             oo.close()

#             oo = open(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/geno_beds_and_intervals/{inds[j]}_batch_{batch}_ploidy_{pl}.intervals", "w")
#             oo.write("\n".join(chromosomes[jj] + ":" + str(starts[jj]) + "-" + str(ends[jj]) for jj in range(len(chromosomes))))
#             oo.close()
             
#### make metadata files for genDB
for batch in sorted(list(set(regions.batch))):
    iis                 = [index for index in range(len(batches)) if batches[index] == batch]
    set_of_ploidy_pairs = pd.DataFrame({"female_ploidies":[list(regions.female_ploidy)[ii] for ii in iis], "male_ploidies":[list(regions.male_ploidy)[ii] for ii in iis]}).drop_duplicates()
    if "M" not in sexes:
        set_of_ploidy_pairs.drop(set_of_ploidy_pairs[set_of_ploidy_pairs['female_ploidies'] == 0].index, inplace = True)
    female_ploidies     = list(set_of_ploidy_pairs.female_ploidies)
    male_ploidies       = list(set_of_ploidy_pairs.male_ploidies)

    for j in range(set_of_ploidy_pairs.shape[0]):
        chromosomes = [chrom_list[ii] for ii in iis if regions.female_ploidy.iloc[ii] == female_ploidies[j] and regions.male_ploidy.iloc[ii] == male_ploidies[j]]
        starts      = [start_list[ii] for ii in iis if regions.female_ploidy.iloc[ii] == female_ploidies[j] and regions.male_ploidy.iloc[ii] == male_ploidies[j]]
        ends        = [end_list[ii] for ii in iis if regions.female_ploidy.iloc[ii] == female_ploidies[j] and regions.male_ploidy.iloc[ii] == male_ploidies[j]]

        oo = open(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/geno_beds_and_intervals/genDB_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}.intervals", "w")
        oo.write("\n".join(chromosomes[jj] + ":" + str(starts[jj]) + "-" + str(ends[jj]) for jj in range(len(chromosomes))))
        oo.close()

        #### make metadata files for genotyping
        batch_lengths     = [end_list[ii] - start_list[ii] + 1 for ii in iis if regions.female_ploidy.iloc[ii] == female_ploidies[j] and regions.male_ploidy.iloc[ii] == male_ploidies[j]] # length of contigs within the batch and ploidy category
        batch_chromosomes = [chrom_list[ii] for ii in iis if regions.female_ploidy.iloc[ii] == female_ploidies[j] and regions.male_ploidy.iloc[ii] == male_ploidies[j]]
        batch_regions     = [regions_list[ii] for ii in iis if regions.female_ploidy.iloc[ii] == female_ploidies[j] and regions.male_ploidy.iloc[ii] == male_ploidies[j]]
        batch_starts      = [start_list[ii] for ii in iis if regions.female_ploidy.iloc[ii] == female_ploidies[j] and regions.male_ploidy.iloc[ii] == male_ploidies[j]]
        batch_ends        = [end_list[ii] for ii in iis if regions.female_ploidy.iloc[ii] == female_ploidies[j] and regions.male_ploidy.iloc[ii] == male_ploidies[j]]

        prev_done_files      = []
        files_to_concatenate = []
        short_segments       = []
        short_segment_starts = []
        short_segment_ends   = []

        chrs = []; sts = []; ens = []
        for k in range(len(batch_lengths)): ## if length of a region longer than length_of_chunk, then cut region into chunks; add smaller regions (if they exist) into short_segment files and genotype together
            if batch_lengths[k] > length_of_chunk:
                tmp_sts = [batch_starts[k]+length_of_chunk*jj for jj in list(range(batch_lengths[k]//length_of_chunk+1))]
                tmp_ens = [batch_starts[k]+length_of_chunk*jj - 1 for jj in list(range(1, batch_lengths[k]//length_of_chunk+1))] + [batch_ends[k]]
                for h in range(len(tmp_sts)):
                    chrs.append(batch_chromosomes[k])
                    sts.append(tmp_sts[h])
                    ens.append(tmp_ens[h])

            else:
                short_segments.append(batch_chromosomes[k])
                short_segment_starts.append(batch_starts[k])
                short_segment_ends.append(batch_ends[k])

        if len(chrs)>0:
            oo = open(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/geno_beds_and_intervals/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_subbatches.intervals", "w")
            oo.write("\n".join(chrs[sb] + ":" + str(sts[sb]) + "-" + str(ens[sb]) for sb in range(len(sts))))
            oo.close()
        
        if len(short_segments)>0:
            oo = open(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/geno_beds_and_intervals/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments.intervals", "w")
            oo.write("\n".join(short_segments[sb] + ":" + str(short_segment_starts[sb]) + "-" + str(short_segment_ends[sb]) for sb in range(len(short_segment_starts))))
            oo.close()