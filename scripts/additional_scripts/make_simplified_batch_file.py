import os
import sys
import pandas as pd

#### load data
regions = pd.read_table(sys.argv[1])
regions = regions.sort_values(by=["batch", "chrom", "start"])
ref_folder = sys.argv[2]

batches = list(regions.batch)
female_ploidy = list(regions.female_ploidy)
male_ploidy = list(regions.male_ploidy)


#### make batches metadata file if it doesn't exist
lines = []
for i in range(regions.shape[0]):
    if str(batches[i])+"\t"+str(female_ploidy[i])+"\t"+str(male_ploidy[i])+"\n" not in lines:
        lines.append(str(batches[i])+"\t"+str(female_ploidy[i])+"\t"+str(male_ploidy[i])+"\n")

oo = open("~/"+ref_folder+"/ref/batches_"+ref_folder+".txt", "w")
oo.write("batch\tfemale_ploidy\tmale_ploidy\n")
for line in lines:
    oo.write(line)
oo.close()
