import os
import numpy as np
import subprocess
from gwf import Workflow, AnonymousTarget
from templates import *
import pandas as pd


stragglers = []
gwf = Workflow()

info = pd.read_csv('primate_ref_samples.csv', sep=',')
species={}
for _, row in info.iterrows():
    if row.FOLDER not in species:
        species[row.FOLDER]=row.REFERENCE_FOLDER
    else:
        continue


for specie in species.keys():
    if os.path.isfile("/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/filteredVCF/bcf_step1/{specie}_all_chr.sorted.bcf".format(specie=specie)):
        # Merge mask
        jobid_merge_mask = "merge_mask_" + specie
        gwf.target_from_template(jobid_merge_mask,
                            merge_mask(path="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/filteredVCF/pos_bed_cov_based/{specie}_batch".format(specie=specie),
                                    out='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/filteredVCF/pos_bed_cov_based/{specie}_all_chr.bed'.format(specie=specie),
                                    done="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/done/{jobid_merge_mask}".format(specie=specie, jobid_merge_mask=jobid_merge_mask)))
        
        # Remove contigs <1mb long
        jobid_remove_small = 'remove_small_contigs_' + specie
        gwf.target_from_template(jobid_remove_small,
                                remove_small_contigs(
                                    bcf="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/filteredVCF/bcf_step1/{specie}_all_chr.sorted.bcf".format(specie=specie),
                                    bed='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/contigs_longer_than1000000.bed'.format(specie=specie),
                                    out='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/filteredVCF/bcf_step1/{specie}_all_chr.1mb.sorted.bcf'.format(specie=specie),
                                    done="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/done/{jobid_remove_small}".format(specie=specie, jobid_remove_small=jobid_remove_small)
                                )
        )


        # Detect the runs of homozygosity 
        jobid_run_roh = "ROH_" + specie
        gwf.target_from_template(jobid_run_roh,
                            run_roh(sorted_bcf="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/filteredVCF/bcf_step1/{specie}_all_chr.1mb.sorted.bcf".format(specie=specie),
                                    out_dir='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/stats/ROH/{specie}'.format(specie=specie),
                                    hom_file='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/stats/ROH/{specie}.hom'.format(specie=specie),
                                    done="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/done/{jobid_run_roh}".format(specie=specie, jobid_run_roh=jobid_run_roh)))
        

    # Handle ROH file, add col of callability, calculate froh and number of windows >= 1mb 
        jobid_get_roh_stats = "ROH_stats_" + specie
        gwf.target_from_template(jobid_get_roh_stats,  # Fixed: was using jobid_run_roh instead
                            get_roh_stats(roh_file='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/stats/ROH/{specie}.hom'.format(specie=specie),
                                        region_file= '/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref}/ref/regions_{ref}_updated.txt'.format(ref=species[specie]),
                                        sex_file = '/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/{species}_individuals.txt'.format(species=specie.split('_')[0]),
                                        species_ssp='{specie}'.format(specie=specie),
                                        out='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/stats/ROH/{specie}_roh_stats.csv'.format(specie=specie), 
                                        figure_seg='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/stats/ROH/{specie}_seg_figure.pdf'.format(specie=specie),
                                        figure_dist='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/stats/ROH/{specie}_dist_figure.pdf'.format(specie=specie),
                                        mask='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/filteredVCF/pos_bed_cov_based/{specie}_all_chr.bed'.format(specie=specie),
                                        done='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/done/{jobid_get_roh_stats}'.format(specie=specie, jobid_get_roh_stats=jobid_get_roh_stats),
                                        done_prev="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/done/{jobid_run_roh}".format(specie=specie, jobid_run_roh=jobid_run_roh)))
        
        # get Het stats
        jobid_het_stats = 'Het_Stats_'+ specie
        gwf.target_from_template(jobid_het_stats,
                                get_het_stats(
                                    bcf_file="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/filteredVCF/bcf_step1/{specie}_all_chr.1mb.sorted.bcf".format(specie=specie),
                                    ploidy_file='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{ref}/ref/regions_{ref}_updated.txt'.format(ref=species[specie]),
                                    mask_file='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/filteredVCF/pos_bed_cov_based/{specie}_all_chr.bed'.format(specie=specie),
                                    sexes_file='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/{species}_individuals.txt'.format(species=specie.split('_')[0]),
                                    out='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/stats/HET/{specie}_het_stats.csv'.format(specie=specie),
                                    done='/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{specie}/done/{jobid_het_stats}'.format(specie=specie, jobid_het_stats=jobid_het_stats),
                                )

        )
    else:
        stragglers.append(specie)

print('species not processed: ', stragglers)
# cat ../../../data/gVCFs_recalling_10_12_2024/Gorilla_beringei_ssp/filteredVCF/pos_bed_cov_based/Gorilla_beringei_ssp_batch* | sort -k 1,1 -k2,2n | bedtools merge > ../../../data/gVCFs_recalling_10_12_2024/Gorilla_beringei_ssp/filteredVCF/pos_bed_cov_based/Gorilla_beringei_ssp_all_chr.bed
