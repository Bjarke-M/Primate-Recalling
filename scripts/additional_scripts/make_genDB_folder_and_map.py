import os
import sys
import pandas as pd
import subprocess

#### load data
group = sys.argv[1]
batch = sys.argv[2]
inds = str(sys.argv[3]).split(",")
ploidies = str(sys.argv[4]).split(",")
folder_name = sys.argv[5]

#### create folders
subprocess.run(f'''mkdir -p {folder_name}/''', shell=True, text=True)

#### create cohort file
pd.DataFrame({"IND_ID" : inds,
              "path"   : [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/{inds[ii]}_batch_{batch}_ploidy_{ploidies[ii]}.gvcf.gz" for ii in range(len(inds))]}).to_csv(f"{folder_name}/cohort.sample_map", sep='\t', header=False, index=False)
