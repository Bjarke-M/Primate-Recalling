import os
import sys
import time

#### load data
group = str(sys.argv[1])
job_name = str(sys.argv[2])
path_to_run_accession = str(sys.argv[3])
expected_md5  = str(sys.argv[4])

#### run
os.system(f"md5sum {path_to_run_accession} > {path_to_run_accession}.md5")

#### get observed md5 value and compare
file = open(f"{path_to_run_accession}.md5", "r")
ln = file.readline()
observed_md5 = ln.split("/")[0].rstrip()

if expected_md5 == observed_md5:
    os.system(f"""
            echo "Download succesful: {path_to_run_accession.split("/")[-1]}" > {path_to_run_accession}.message

            touch ~/{group}/done/{job_name}
            """)
else:
    os.system(f"""
            echo "Download failed: {path_to_run_accession.split("/")[-1]}" > {path_to_run_accession}.message
            """)

