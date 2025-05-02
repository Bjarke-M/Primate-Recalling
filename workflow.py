"""
------------------------------------------------------------------------------------------------------------------------
This workflow handles raw files (reference genomes, paired-end fastqs), maps reads, calls and genotypes variants.
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Authors: Juraj Bergman, Vasili Pankratov, Bjarke M. Pedersen
Date: 26/03/2025
------------------------------------------------------------------------------------------------------------------------
"""

import os
import numpy as np
import subprocess
from gwf import Workflow, AnonymousTarget
from templates import *
import pandas as pd
gwf = Workflow()

## load global metadata files
user       = "bjarkemp"
genus      = "Papio"
# data       = pd.read_table(f"/faststorage/project/primatediversity/people/{user}/workflow_recalling/{genus}_with_Kuderna_IDs.txt")
# references = pd.read_table(f"/faststorage/project/primatediversity/people/{user}/workflow_recalling/{genus}_references.txt")
data_path = f"/faststorage/project/primatediversity/people/{user}/workflow_recalling/{genus}.txt"
data       = pd.read_table(f"/faststorage/project/primatediversity/people/{user}/workflow_recalling/{genus}.txt")
references = pd.read_table(f"/faststorage/project/primatediversity/people/{user}/workflow_recalling/{genus}_references.txt")

## make data frame that contains names of species-specific folders and the reference folders used to map the species
species_and_refs = pd.DataFrame({"FOLDER": data.FOLDER, "REFERENCE_FOLDER": data.REFERENCE_FOLDER}).drop_duplicates()
species_and_refs = species_and_refs.reset_index(drop=True)

## merge dataframes to have species-specific folder, reference folder and fastq ftps in same dataframe
species_and_refs = species_and_refs.merge(references, how = "left")

"""
------------------------------------------------------------------------------------------------------------------------
A. REFERENCE-ASSOCIATED JOBS
------------------------------------------------------------------------------------------------------------------------
"""

## run jobs that download reference, make a reference that contains only contigs more or equal to 1000 bp in length and makes a referece-specific regions file (to be manually modified as necessary due to X, PAR, Y, mitochondrial ploidy, etc.)

## submit jobs
# for i in range(references.shape[0]):
#     ## A.0. Initialize directories
#     subprocess.run(["mkdir", "-p", "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done".format(references.REFERENCE_FOLDER[i])])
#     subprocess.run(["mkdir", "-p", "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref".format(references.REFERENCE_FOLDER[i])])
#     subprocess.run(["mkdir", "-p", "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/fastq".format(references.REFERENCE_FOLDER[i])])
#     subprocess.run(["mkdir", "-p", "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/bam".format(references.REFERENCE_FOLDER[i])])
#     subprocess.run(["mkdir", "-p", "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/gVCF".format(references.REFERENCE_FOLDER[i])])

#     ## A.1. download reference
#     jobid_download_ref = "download_ref_" + references.ref_genome_name[i].replace("-","_")
#     gwf.target_from_template(jobid_download_ref,
#                              download_ref(ftp  = references.ftp[i],
#                                           out  = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-3]),
#                                           done = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i],jobid_download_ref)))

#     ## A.2. mask reference - uncomment if reference needs to be masked; "masked_regions.bed is expected in the reference folder
#     # jobid_mask_reference = "mask_reference_" + references.ref_genome_name[i].replace("-", "_")
#     # gwf.target_from_template(jobid_mask_reference,
#     #                          mask_reference(input="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-3]),
#     #                                         bed="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref/{}".format(references.REFERENCE_FOLDER[i], "masked_regions.bed"),
#     #                                         output="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-7] + "_masked.fna"),
#     #                                         done_prev="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_download_ref),
#     #                                         done="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_mask_reference)))

#     ## A.3. cut out contigs less than 1000 bp in length
#     jobid_cut_contigs = "cut_contigs_" + references.ref_genome_name[i].replace("-", "_")
#     gwf.target_from_template(jobid_cut_contigs,
#                              cut_contigs(input     = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-3]),
#                              # cut_contigs(input      = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-7]+"_masked.fna"), # if masked reference present comment out previous line and uncomment this line
#                                          output     = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-7]+"_LargerThan1000bp.fasta"),
#                                          min_length = 1000,
#                                          done_prev  = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_download_ref),
#                                          # done_prev="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_mask_reference), # if masked reference present comment out previous line and uncomment this line
#                                          done       = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_cut_contigs)))

#     ## A.4. make reference-associated files
#     jobid_make_fasta = "make_fasta_" + references.ref_genome_name[i].replace("-","_")
#     gwf.target_from_template(jobid_make_fasta,
#                              make_fasta(input     = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-3]),
#                              # make_fasta(input="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/ref/{}".format(references.REFERENCE_FOLDER[i], str(references.ftp[i]).split("/")[-1][:-7]+"_masked.fna"), ## if masked reference present comment out previous line and uncomment this line
#                                         done_prev = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_cut_contigs),
#                                         done      = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_make_fasta)))

#     ## A.5. make regions files
#     jobid_make_regions = "make_regions_" + references.ref_genome_name[i].replace("-","_")
#     gwf.target_from_template(jobid_make_regions,
#                              make_regions(refFolder = references.REFERENCE_FOLDER[i],
#                                           input     = str(references.ftp[i]).split("/")[-1][:-3],
#                                           # input=str(references.ftp[i]).split("/")[-1][:-7]+"_masked.fna", ## if masked reference present comment out previous line and uncomment this line
#                                           done_prev = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_make_fasta),
#                                           done      = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{}/done/{}".format(references.REFERENCE_FOLDER[i], jobid_make_regions)))

"""
------------------------------------------------------------------------------------------------------------------------
B. SAMPLE DOWNLOAD AND PREPARATION FOR MAPPING
------------------------------------------------------------------------------------------------------------------------
"""
## submit jobs
# for i in range(species_and_refs.shape[0]):
#     group = species_and_refs.FOLDER[i]
#     #if group == "Macaca_mulatta_ssp":
#     inds  = list(data.loc[data.FOLDER == group].IND_ID.drop_duplicates())
#     srrs  = list(data.loc[data.FOLDER == group].fastq_ftp)

#     ## B.0. Initialize directories
#     subprocess.run(["mkdir", "-p", f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done"])
#     subprocess.run(["mkdir", "-p", f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq"])
#     subprocess.run(["mkdir", "-p", f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam"])
#     subprocess.run(["mkdir", "-p", f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF"])

#     for ind in inds:
#         ## B.1. Download fastq files per individual 
#         data_subset = data.loc[data.IND_ID == ind]
#         jobid_download_per_individual = f"download_per_individual_{group}_{ind}"
#         gwf.target_from_template(jobid_download_per_individual,
#                                 download_per_individual(group          = group,
#                                                         ind            = ind,
#                                                         run_accessions = ",".join(list(data_subset.fastq_ftp)),
#                                                         md5s           = ",".join(list(data_subset.fastq_md5)),
#                                                         done           = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_download_per_individual}"))

#         ## B.2. Concatenate fastq files per individual if more than one accession are available, else rename the fastq files; prepare uBAMs for mapping
#         if data_subset.shape[0] > 2:
#             jobid_concat_or_rename_fastqs = f"concat_or_rename_fastqs_{group}_{ind}"
#             gwf.target_from_template(jobid_concat_or_rename_fastqs,
#                                     concatfastqs(group     = group,
#                                                 ind       = ind,
#                                                 fastqs_1  = ["/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/" + group + "/fastq/" + data_subset.loc[data_subset.R1_or_R2 == "R1"].fastq_ftp.iloc[jj].split("/")[-1] for jj in range(data_subset.loc[data_subset.R1_or_R2 == "R1"].shape[0])],
#                                                 fastqs_2  = ["/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/" + group + "/fastq/" + data_subset.loc[data_subset.R1_or_R2 == "R2"].fastq_ftp.iloc[jj].split("/")[-1] for jj in range(data_subset.loc[data_subset.R1_or_R2 == "R2"].shape[0])],
#                                                 prev_done = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_download_per_individual}",# +["/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/" + group + "/done/download_pe2_" + group + "_" + srr.replace("/", "_") for srr in data_subset.fastq_ftp],
#                                                 done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_concat_or_rename_fastqs}"))
#         else:
#             jobid_concat_or_rename_fastqs = f"concat_or_rename_fastqs_{group}_{ind}"
#             gwf.target_from_template(jobid_concat_or_rename_fastqs,
#                                     renamefastqs(group     = group,
#                                                 ind       = ind,
#                                                 fastq_1   = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/" + group + "/fastq/" + data_subset.loc[data_subset.R1_or_R2 == "R1"].fastq_ftp.iloc[0].split("/")[-1],
#                                                 fastq_2   = "/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/" + group + "/fastq/" + data_subset.loc[data_subset.R1_or_R2 == "R2"].fastq_ftp.iloc[0].split("/")[-1],
#                                                 prev_done = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_download_per_individual}",# +["/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/" + group + "/done/download_pe2_" + group + "_" + srr.replace("/", "_") for srr in data_subset.fastq_ftp],
#                                                 done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_concat_or_rename_fastqs}"))

#         ## B.3. Make uBAMs from concatenated or renamed fastqs
#         jobid_makeuBAM = f"makeuBAM_{group}_{ind}"
#         gwf.target_from_template(jobid_makeuBAM,
#                                 makeuBAM(group=group,
#                                         ind       = ind,
#                                         fastq1    = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq/{ind}_R1.fastq.gz",
#                                         fastq2    = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq/{ind}_R2.fastq.gz",
#                                         prev_done = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_concat_or_rename_fastqs}",
#                                         done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_makeuBAM}"))

#         ## B.4. Split uBAMs
#         jobid_splituBAM = f"splituBAM_{group}_{ind}"
#         gwf.target_from_template(jobid_splituBAM,
#                                 splituBAM(group     = group,
#                                         ind       = ind,
#                                         prev_done = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_makeuBAM}"],
#                                         done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_splituBAM}"))

"""
------------------------------------------------------------------------------------------------------------------------
C. MAPPING
------------------------------------------------------------------------------------------------------------------------
"""
## prequisite to start mapping is the existence of the f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_splituBAM}" done file

## submit jobs
# for i in range(species_and_refs.shape[0]):
#     group      = species_and_refs.FOLDER[i]
#     ref_folder = species_and_refs.REFERENCE_FOLDER[i]
#     ref_path   = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/" + [file for file in os.listdir(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/") if file.endswith("_LargerThan1000bp.fasta")][0]
#     inds       = list(data.loc[data.FOLDER == group].IND_ID.drop_duplicates())

#     for ind in inds:
#         jobid_splituBAM = f"splituBAM_{group}_{ind}" ## prequisite to start mapping is the existence of the f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_splituBAM}" done file
#         if os.path.isfile(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_splituBAM}"): ## checks if f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_splituBAM}" done file exists
#             with open(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_nsplitubams.txt") as f:
#                 for l in f:
#                     nbams = int(l.strip())

#             for ishard in range(nbams):
#                 shard = shardstr(ishard)

#                 ## C.1. Mark the adapters in each shard uBAM file
#                 jobid_markadapt = f"markadapt_{group}_{ind}_{shard}"
#                 gwf.target_from_template(jobid_markadapt,
#                                             markadapt(group   = group,
#                                                     ind       = ind,
#                                                     shard     = shard,
#                                                     prev_done = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_splituBAM}"],
#                                                     done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_markadapt}"))

#                 ## C.2. Map each shard uBAM file
#                 jobid_mapBAM = f"mapBAM_{group}_{ind}_{shard}"
#                 gwf.target_from_template(jobid_mapBAM,
#                                             mapBAM(group  = group,
#                                                 ind       = ind,
#                                                 shard     = shard,
#                                                 ref       = ref_path,
#                                                 prev_done = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_markadapt}"],
#                                                 done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_mapBAM}"))

#             ## C.3 Merge mapped shard bams into single bam per individual
#             jobid_mergeBAMs = f"mergeBAMs_{group}_{ind}"
#             gwf.target_from_template(jobid_mergeBAMs,
#                                         mergeBAMs(group     = group,
#                                                 ind       = ind,
#                                                 bams      = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}/shard_{shardstr(ishard)}_markadapt_mapped.bam" for ishard in range(nbams)],
#                                                 prev_done = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/mapBAM_{group}_{ind}_{shardstr(ishard)}" for ishard in range(nbams)],
#                                                 done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_mergeBAMs}"))

#             ## C.4. Mark and remove duplicates from bam
#             jobid_markduplicates = f"markduplicates_{group}_{ind}"
#             gwf.target_from_template(jobid_markduplicates,
#                                         markduplicates(group  = group,
#                                                     ind       = ind,
#                                                     prev_done = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_mergeBAMs}"],
#                                                     done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_markduplicates}"))

#             ## C.5. Sort bam by coordinates
#             jobid_coordsort = f"coordsort_{group}_{ind}"
#             gwf.target_from_template(jobid_coordsort,
#                                         coordsort(group   = group,
#                                                 ind       = ind,
#                                                 prev_done = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_markduplicates}"],
#                                                 done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_coordsort}"))

#             ## C.6. Get coverage statistics
#             regions      = pd.read_table(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/regions_{ref_folder}.txt")
#             regions_list = list(regions.region)
#             chrom_list   = list(regions.chrom)
#             start_list   = [jj + 1 for jj in list(regions.start)]
#             end_list     = list(regions.end)

#             jobid_cov = f"cov_{group}_{ind}"
#             gwf.target_from_template(jobid_cov,
#                                         cov(group       = group,
#                                             ind         = ind,
#                                             regions     = regions_list,
#                                             chromosomes = chrom_list,
#                                             starts      = start_list,
#                                             ends        = end_list,
#                                             prev_done   = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_coordsort}"],
#                                             done        = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_cov}"))


"""
------------------------------------------------------------------------------------------------------------------------
D. CALLING AND GENOTYPING
------------------------------------------------------------------------------------------------------------------------
"""
## pre-calling jobs
for i in range(references.shape[0]):
    ## D.1. find X contigs and update regions file 
    jobid_find_chrX = f'find_chrX_{references.REFERENCE_FOLDER.iloc[i]}'
    gwf.target_from_template(jobid_find_chrX,
                                find_chrX(subset_file = data_path,
                                        ref_folder   = references.REFERENCE_FOLDER.iloc[i],
                                        contigs      = "/faststorage/project/primatediversity/people/vasili/workflow_recalling/special_contigs",
                                        minlen       = 1e6,
                                        prev_done    = [],
                                        done         = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{references.REFERENCE_FOLDER.iloc[i]}/done/{jobid_find_chrX}"))

    ## D.2. make simplified batch file from regions file
    jobid_make_simplified_batch_file = f'make_simplified_batch_file_{references.REFERENCE_FOLDER.iloc[i]}'
    gwf.target_from_template(jobid_make_simplified_batch_file,
                                make_simplified_batch_file(regions_file = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{references.REFERENCE_FOLDER.iloc[i]}/ref/regions_{references.REFERENCE_FOLDER.iloc[i]}_updated.txt",
                                                            ref_folder   = references.REFERENCE_FOLDER.iloc[i],
                                                            prev_done    = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{references.REFERENCE_FOLDER.iloc[i]}/done/{jobid_find_chrX}"],
                                                            done         = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{references.REFERENCE_FOLDER.iloc[i]}/done/{jobid_make_simplified_batch_file}"))

## prequisite to start calling and genotyping is the existence of the f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_cov}" and f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_loc_metadata}" done files
for i in range(species_and_refs.shape[0]):
    group          = species_and_refs.FOLDER[i]
    ref_folder      = species_and_refs.REFERENCE_FOLDER[i]
    ref_path        = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/" + [file for file in os.listdir(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/") if file.endswith("_LargerThan1000bp.fasta")][0]
    inds_to_include = list(data.loc[data.FOLDER == group]["IND_ID"].drop_duplicates())
    inds_and_sexes  = pd.read_table(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/samples_coverage_stats.txt")
    inds_and_sexes  = inds_and_sexes[inds_and_sexes.IND_ID.isin(inds_to_include)]
    inds            = list(inds_and_sexes.IND_ID)
    sexes           = list(inds_and_sexes.gSEX)

    regions = pd.read_table(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/batches_{ref_folder}.txt")
    batches = list(regions.batch)

    ## D.3. make metadata files
    length_of_chunk = 2000000
    jobid_make_simplified_batch_file = f'make_simplified_batch_file_{species_and_refs.REFERENCE_FOLDER.iloc[i]}'
    jobid_make_loc_metadata = f'make_loc_metadata_{group}'
    gwf.target_from_template(jobid_make_loc_metadata,
                                make_batch_metadata(loc          = length_of_chunk,
                                                    group        = group,
                                                    regions_file = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/regions_{ref_folder}_updated.txt",
                                                    inds         = ",".join(inds),
                                                    sexes        = ",".join(sexes),
                                                    ref_folder   = ref_folder,
                                                    prev_done    = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{species_and_refs.REFERENCE_FOLDER.iloc[i]}/done/{jobid_make_simplified_batch_file}"],
                                                    done         = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_loc_metadata}"))

#     ## D.4. Call variants per individual considering region ploidy (as determined by the sex of the individual)
    for j in range(len(inds)):
        jobid_cov = f"cov_{group}_{inds[j]}" ## prequisite to start calling and genotyping is the existence of the f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_cov}" done file
        if os.path.isfile(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_cov}"): ## checks if  f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_cov}" done file exists
            ## get ploidy of regions based on sex of individual
            if sexes[j] == "M":
                ploidies_list = list(regions.male_ploidy)
            else:
                ploidies_list = list(regions.female_ploidy)

            for batch in sorted(list(set(regions.batch))):
                iis = [index for index in range(len(batches)) if batches[index] == batch]
                set_of_ploidies = sorted(list(set([ploidies_list[ii] for ii in iis if ploidies_list[ii] != 0])))
                for pl in set_of_ploidies:
                    jobid_call = f"call_with_bed_{group}_{inds[j]}_batch_{batch}_ploidy_{pl}"
                    gwf.target_from_template(jobid_call,
                                                call_batch_with_bed(group = group,
                                                    ind                = inds[j],
                                                    batch              = batch,
                                                    bed                = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/beds_and_intervals/{inds[j]}_batch_{batch}_ploidy_{pl}.bed",
                                                    intervals          = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/beds_and_intervals/{inds[j]}_batch_{batch}_ploidy_{pl}.intervals",
                                                    ref                = ref_path,
                                                    ploidy             = pl,
                                                    prev_done          = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_cov}", f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_loc_metadata}"],
                                                    done               = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_call}"))

#     ## D.5. Make a GenomicsDB folder and cohort.sample_map file for each batch and ploidy category within batch
    for batch in sorted(list(set(regions.batch))):
        iis                 = [index for index in range(len(batches)) if batches[index] == batch]
        set_of_ploidy_pairs = pd.DataFrame({"female_ploidies":[list(regions.female_ploidy)[ii] for ii in iis], "male_ploidies":[list(regions.male_ploidy)[ii] for ii in iis]}).drop_duplicates()
        if "M" not in sexes:
            set_of_ploidy_pairs.drop(set_of_ploidy_pairs[set_of_ploidy_pairs['female_ploidies'] == 0].index, inplace = True)
        female_ploidies     = list(set_of_ploidy_pairs.female_ploidies)
        male_ploidies       = list(set_of_ploidy_pairs.male_ploidies)

        for j in range(set_of_ploidy_pairs.shape[0]):
            if female_ploidies[j] == male_ploidies[j]:
                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{female_ploidies[j]}'
                gwf.target_from_template(jobid_make_genDB_folder_and_map,
                                            make_genDB_folder_and_map(group       = group,
                                                                    batch       = batch,
                                                                    inds        = ",".join(inds),
                                                                    ploidies    = ",".join(len(inds)*[str(female_ploidies[j])]),
                                                                    folder_name = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}/",
                                                                    prev_done   = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/call_with_bed_{group}_{inds[jj]}_batch_{batch}_ploidy_{female_ploidies[j]}" for jj in range(len(inds))],
                                                                    done        = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_genDB_folder_and_map}"))

            elif female_ploidies[j] == 2 and male_ploidies[j] == 1:
                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{group}_batch_{batch}_fploidy_2_mploidy_1'
                gwf.target_from_template(jobid_make_genDB_folder_and_map,
                                            make_genDB_folder_and_map(group       = group,
                                                                    batch       = batch,
                                                                    inds        = ",".join(inds),
                                                                    ploidies    = ",".join(["2" if sexes[jj] == "F" else "1" for jj in range(len(sexes))]),
                                                                    folder_name = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}/",
                                                                    prev_done   = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/call_with_bed_{group}_{inds[jj]}_batch_{batch}_ploidy_2" if sexes[jj] == "F" else f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/call_with_bed_{group}_{inds[jj]}_batch_{batch}_ploidy_1" for jj in range(len(inds))],
                                                                    done        = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_genDB_folder_and_map}"))

            elif female_ploidies[j] == 0 and male_ploidies[j] == 1:
                jobid_make_genDB_folder_and_map = f'make_genDB_folder_and_map_{group}_batch_{batch}_fploidy_0_mploidy_1'
                gwf.target_from_template(jobid_make_genDB_folder_and_map,
                                            make_genDB_folder_and_map(group       = group,
                                                                    batch       = batch,
                                                                    inds        = ",".join([inds[jj] for jj in range(len(inds)) if sexes[jj] == "M"]),
                                                                    ploidies    = ",".join(sexes.count("M")*["1"]),
                                                                    folder_name = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_0_mploidy_1/",
                                                                    prev_done   = [f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/call_with_bed_{group}_{inds[jj]}_batch_{batch}_ploidy_1" for jj in range(len(inds)) if sexes[jj] == "M"],
                                                                    done        = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_genDB_folder_and_map}"))

            ## D.6. Make a GenomicsDB workspace for each batch and ploidy category within batch
            jobid_make_genDB = f'make_genDB_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
            gwf.target_from_template(jobid_make_genDB,
                                        make_genDB_with_bed(group       = group,
                                                            batch       = batch,
                                                            fploidy     = female_ploidies[j],
                                                            mploidy     = male_ploidies[j],
                                                            intervals   = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/beds_and_intervals/genDB_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}.intervals",
                                                            prev_done   = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_genDB_folder_and_map}",
                                                            done        = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_genDB}"))

            ## D.7. Genotype for each batch and ploidy category within batch; separate regions into chunks of length_of_chunk if the are longer than length_of_chunk               
            prev_done_files = []
            files_to_concatenate = []
            if os.path.isfile(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/beds_and_intervals/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_subbatches.intervals"):
                subbatch_file = pd.read_table(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/beds_and_intervals/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_subbatches.intervals", header = None)
                for sb in range(subbatch_file.shape[0]):
                    ## Genotype
                    jobid_GenotypeGVCFs = f'GenotypeGVCFs_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_subbatch_{sb}'
                    prev_done_files.append(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_GenotypeGVCFs}")
                    files_to_concatenate.append(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_subbatch_{sb}_gt.gvcf.gz")
                    if False == os.path.isfile(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/IndexGVCFs_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}"): ## checks if the Index done file exists (if it does exist it is an indication that the genotyping was completed for this batch and and ploidy category)
                        gwf.target_from_template(jobid_GenotypeGVCFs,
                                                GenotypeGVCFs_subbatch_with_bed(group      = group,
                                                                                batch      = batch,
                                                                                fploidy    = female_ploidies[j],
                                                                                mploidy    = male_ploidies[j],
                                                                                subbatch   = str(sb),
                                                                                interval   = subbatch_file.iloc[sb,0],
                                                                                out        = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_subbatch_{sb}_gt.gvcf.gz",
                                                                                ref        = ref_path,
                                                                                cores      = 1,
                                                                                memory     = "32g",
                                                                                prev_done  = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_genDB}",
                                                                                done       = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_GenotypeGVCFs}"))

            if os.path.isfile(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/beds_and_intervals/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments.intervals"): ## genotype smaller regions shorter than length_of_chunk (if they exist)
                ## Genotype
                jobid_GenotypeGVCFs = f'GenotypeGVCFs_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_short_segments'
                prev_done_files.append(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_GenotypeGVCFs}")
                files_to_concatenate.append(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_short_segments_gt.gvcf.gz")
                if False == os.path.isfile(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/IndexGVCFs_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}"): ## checks if the Index done file exists (if it does exist it is an indication that the genotyping was completed for this batch and and ploidy category)
                    gwf.target_from_template(jobid_GenotypeGVCFs,
                                            GenotypeGVCFs_with_bed(group         = group,
                                                                    batch       = batch,
                                                                    fploidy     = female_ploidies[j],
                                                                    mploidy     = male_ploidies[j],
                                                                    intervals   = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/beds_and_intervals/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_loc_{length_of_chunk}_short_segments.intervals",
                                                                    out         = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_short_segments_gt.gvcf.gz",
                                                                    ref         = ref_path,
                                                                    cores       = 1,
                                                                    memory      = "8g",
                                                                    prev_done   = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_make_genDB}",
                                                                    done        = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_GenotypeGVCFs}"))

            if len(files_to_concatenate) > 1: ## if chunking was done, concatenate all chunks within a batch and ploidy category, then index the concatenated gvcf and finally remove individual-based GVCFs
                ## Concatenate
                jobid_picardconcat = f'picardconcat_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
                gwf.target_from_template(jobid_picardconcat,
                                            picardconcat(vcfs      = files_to_concatenate,
                                                        vcf       = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_gt.gvcf.gz",
                                                        prev_done = prev_done_files,
                                                        done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_picardconcat}"))
                
                ## Index
                jobid_index = f'IndexGVCFs_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
                gwf.target_from_template(jobid_index,
                                        IndexGVCFs(group      = group,
                                                    batch     = batch,
                                                    fploidy   = female_ploidies[j],
                                                    mploidy   = male_ploidies[j],
                                                    prev_done = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_picardconcat}",
                                                    done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_index}"))
                
            else:  ## if chunking was NOT done, rename the gvcf to follow convention, then index the renamed gvcf and finally remove individual-based GVCFs
                ## Rename
                jobid_renameGVCF = f'renameGVCF_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
                gwf.target_from_template(jobid_renameGVCF,
                                            renameGVCF(vcf_in    = files_to_concatenate[0],
                                                        vcf_out   = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}_gt.gvcf.gz",
                                                        prev_done = prev_done_files,
                                                        done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_renameGVCF}"))
                
                ## Index
                jobid_index = f'IndexGVCFs_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
                gwf.target_from_template(jobid_index,
                                        IndexGVCFs(group     = group,
                                                    batch     = batch,
                                                    fploidy   = female_ploidies[j],
                                                    mploidy   = male_ploidies[j],
                                                    prev_done = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_renameGVCF}",
                                                    done      = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_index}"))

"""
------------------------------------------------------------------------------------------------------------------------
E. CLEAN UP
------------------------------------------------------------------------------------------------------------------------
"""
# for i in range(species_and_refs.shape[0]):
#     group          = species_and_refs.FOLDER[i]
#     ref_folder     = species_and_refs.REFERENCE_FOLDER[i]
#     ref_path       = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/" + [file for file in os.listdir(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/") if file.endswith("_LargerThan1000bp.fasta")][0]
    
#     inds_to_include = list(data.loc[data.FOLDER == group]["IND_ID"].drop_duplicates())
#     inds_and_sexes  = pd.read_table(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/samples_coverage_stats.txt")
#     inds_and_sexes  = inds_and_sexes[inds_and_sexes.IND_ID.isin(inds_to_include)]
#     inds            = list(inds_and_sexes.IND_ID)
#     sexes           = list(inds_and_sexes.gSEX)

#     regions = pd.read_table(f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref_folder}/ref/regions_{ref_folder}_updated.txt")
#     regions = regions.sort_values(by=["batch", "chrom", "start"])
#     batches = list(regions.batch)

#     ## E.1. Convert bam to cram 
#     for ind in inds:
#         jobid_bam_to_cram = f"bam_to_cram_{group}_{ind}"
#         gwf.target_from_template(jobid_bam_to_cram,
#                                     bam_to_cram(path_to_bam_folder = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam",
#                                             ind                 = ind,
#                                             ref                 = ref_path,
#                                             done                = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_bam_to_cram}"))

#     ## E.2. Remove individual GVCFs
#     for batch in sorted(list(set(regions.batch))):
#         iis                 = [index for index in range(len(batches)) if batches[index] == batch]
#         set_of_ploidy_pairs = pd.DataFrame({"female_ploidies":[list(regions.female_ploidy)[ii] for ii in iis], "male_ploidies":[list(regions.male_ploidy)[ii] for ii in iis]}).drop_duplicates()
#         female_ploidies     = list(set_of_ploidy_pairs.female_ploidies)
#         male_ploidies       = list(set_of_ploidy_pairs.male_ploidies)

#         for j in range(set_of_ploidy_pairs.shape[0]):
#             jobid_index = f'IndexGVCFs_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
#             jobid_removeIndGVCFs = f'removeIndGVCFs_{group}_batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}'
#             gwf.target_from_template(jobid_removeIndGVCFs,
#                                      removeIndGVCFs(cohort_map = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{female_ploidies[j]}_mploidy_{male_ploidies[j]}/cohort.sample_map",
#                                                      prev_done  = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_index}",
#                                                      done       = f"/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/done/{jobid_removeIndGVCFs}"))