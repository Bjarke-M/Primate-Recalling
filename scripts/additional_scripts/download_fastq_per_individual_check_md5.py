import os
import sys
import time

#### load data
group = sys.argv[1]
ind = sys.argv[2]
run_accessions = str(sys.argv[3]).split(",")
md5s = str(sys.argv[4]).split(",")

#### run
attempt = 0
a_files = [f"download_pe2_{group}_" + run_accession.replace("/", "_") for run_accession in run_accessions]
a_exist = [f for f in a_files if os.path.isfile(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/"+f)]
while len(a_files) != len(a_exist):
    os.system(f"""echo Attempt: "{attempt}"
    """)
    for i in range(len(run_accessions)):

        run_accession  = run_accessions[i]
        run_accession_name = run_accessions[i].split("/")[-1]
        run_accession_md5 = md5s[i]
        job_name = f"download_pe2_{group}_" + run_accessions[i].replace("/", "_")

        if job_name in os.listdir(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/"):
            os.system(f"""
            echo "{job_name} present in /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/"
            echo "Skipping..."
            """)
        else:
            os.system(f"""
            echo "Starting download: {run_accession_name}"
            
            rm -fr /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq/{run_accession_name}
            
            wget --progress=dot:giga --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq/ {run_accession}
            
            python /faststorage/project/primatediversity/people/juraj/python_scripts/check_md5.py {group} {job_name} /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq/{run_accession_name} {run_accession_md5}

            cat /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq/{run_accession_name}.message

            rm -f /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq/{run_accession_name}.md5
            rm -f /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq/{run_accession_name}.message

            """)

    a_exist = [f for f in a_files if os.path.isfile(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/"+f)]
    attempt += 1
    l1 = len(a_files)
    l2 = len(a_exist)
    os.system(f"""echo "a_files: {l1}"
        """)
    os.system(f"""echo "a_exist: {l2}"
            """)
    if attempt == 100:
        break

if len(a_files) == len(a_exist):
    os.system(f"""
              echo "Making master done file: /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/download_per_individual_{group}_{ind}"
              touch /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/download_per_individual_{group}_{ind}
              """)
