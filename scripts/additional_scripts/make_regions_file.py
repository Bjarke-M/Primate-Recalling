import os
import sys
import pandas as pd

#### load data
fastaFolder = sys.argv[1] # "Equus_przewalskii"
bed = sys.argv[2] # "GCA_037783145.1_EquPr2_genomic_LargerThan1000bp.chrs"

#### run
chrs = pd.read_table("/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/"+fastaFolder+"/ref/"+bed, header=None)
chr_name = list(chrs.iloc[:,0])
chr_len = list(chrs.iloc[:,1])
temp_group = []; temp_len = 0

chr_groups = []; chr_group_names = []; gn = 0
temp_group = []; temp_len = 0

for ii in range(len(chr_name)):
    if temp_len > 30000000:
        chr_groups.append(temp_group)
        chr_group_names.append(str(gn))
        gn += 1
        temp_group = []; temp_len = 0
        temp_group.append(chr_name[ii])
        temp_len += chr_len[ii]
    else:
        temp_group.append(chr_name[ii])
        temp_len += chr_len[ii]
    if ii == len(chr_name)-1:
        chr_groups.append(temp_group)
        chr_group_names.append(str(gn))

oo = open("/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/"+fastaFolder+"/ref/regions_"+fastaFolder+".txt", "w")
oo.write("region\tchrom\tstart\tend\tbatch\tfemale_ploidy\tmale_ploidy\n")
for ii in range(len(chr_groups)):
    for jj in chr_groups[ii]:
        oo.write(jj + "\t" + jj + "\t0\t"+str(chr_len[chr_name.index(jj)])+ "\t" + str(ii)+"\t2\t2\n")
oo.close()