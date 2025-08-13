"""
------------------------------------------------------------------------------------------------------------------------
This file containes gwf functions of the stats pipeline.
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Authors: Bjarke M. Pedersen, Juraj Bergman, Vasili Pankratov ;)
Date: 29/06/2025
------------------------------------------------------------------------------------------------------------------------
"""

#### Importing and setting default options
from gwf import Workflow, AnonymousTarget
gwf = Workflow()

job_header = '''
    echo "JOBID:" $PBS_JOBID
    echo "NODE :" $HOSTNAME
    echo "USER :" $USER
    source ~/.bashrc
    conda activate /home/bjarkemp/miniforge3/envs/primate_stats
    echo "CONDA:" $CONDA_DEFAULT_ENV
'''

default_options = {"cores" : 1, 'memory': "32g", 'walltime': "02:00:00", 'account': "primatediversity"}

#### Templates
"""
------------------------------------------------------------------------------------------------------------------------
RUNS OF HOMOZYGOSITY ASSOCIATED JOBS
------------------------------------------------------------------------------------------------------------------------
"""

def merge_mask(path,out,done):
    inputs=[]
    outputs=[out,done]
    options = default_options.copy()
    spec = job_header+'''
    cat {path}* | sort -k 1,1 -k2,2n | bedtools merge > {out}
    touch {done}
    '''.format(path=path,out=out,done=done)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def remove_small_contigs(bcf,bed,out,done):
    inputs=[bcf,bed]
    outputs=[out,done]
    options = default_options.copy()
    spec = job_header+'''
    bcftools view --regions-file {bed} {bcf} -Ob -o {out}
    bcftools index {out}
    touch {done}
    '''.format(bcf=bcf,bed=bed,out=out,done=done)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# For BCFtools 

# def run_roh(sorted_bcf, out, done):
#     inputs = [sorted_bcf]  # Include the input BCF file
#     outputs = [out, done]  # Include both output files
#     options = {"cores" : 1, 'memory': "200g", 'walltime': "02:00:00", 'account': "primatediversity"}
    
#     spec = job_header+'''
#     mkdir -p "$(dirname {out})" && touch "{out}"
#     bcftools roh -G30 -o {out} --AF-dflt 0.4 {sorted_bcf}
#     touch {done}
#     '''.format(out=out, sorted_bcf=sorted_bcf, done=done)
    
#     # Use the correct AnonymousTarget constructor
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def get_roh_stats(roh_file, mask, region_file, sex_file, species_ssp,out, figure_seg, figure_dist, tempfile, done, done_prev):
#     inputs = [done_prev, roh_file, mask]  # Also include roh_file as input
#     outputs = [out, figure_seg, figure_dist, done]   # Include all output files
#     options = default_options.copy()
    
#     spec = job_header+'''
#     cat {roh_file} | grep 'RG' > {tempfile}
#     python simple.py {tempfile} {mask} {region_file} {sex_file} {species_ssp} {out} {figure_seg} {figure_dist}
#     rm {tempfile}
#     touch {done}
#     '''.format(roh_file=roh_file, region_file=region_file, sex_file=sex_file, species_ssp=species_ssp,
#                figure_seg=figure_seg, figure_dist=figure_dist, out=out, tempfile=tempfile, done=done, mask=mask)
    
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# Plink Version
def run_roh(sorted_bcf, out_dir, done, hom_file):
    inputs = [sorted_bcf]  # Include the input BCF file
    outputs = [done, hom_file]  # Include both output files
    options = {"cores" : 1, 'memory': "4g", 'walltime': "02:00:00", 'account': "primatediversity"}
    
    spec = job_header+'''
    mkdir -p "$(dirname {out_dir})"
    plink --homozyg --bcf {sorted_bcf} --out {out_dir} --allow-extra-chr
    touch {done}
    '''.format(out_dir=out_dir, sorted_bcf=sorted_bcf, done=done)
    
    # Use the correct AnonymousTarget constructor
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def get_roh_stats(roh_file, mask, region_file, sex_file, species_ssp,out, figure_seg, figure_dist, done, done_prev):
    inputs = [done_prev, roh_file, mask]  # Also include roh_file as input
    outputs = [out, figure_seg, figure_dist, done]   # Include all output files
    options = default_options.copy()
    
    spec = job_header+'''
    python roh_stats.py {roh_file} {mask} {region_file} {sex_file} {species_ssp} {out} {figure_seg} {figure_dist}
    touch {done}
    '''.format(roh_file=roh_file, region_file=region_file, sex_file=sex_file, species_ssp=species_ssp,
               figure_seg=figure_seg, figure_dist=figure_dist, out=out, done=done, mask=mask)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_het_stats(bcf_file, ploidy_file, mask_file, sexes_file, out, done):
    inputs=[bcf_file,ploidy_file,mask_file,sexes_file]
    outputs=[out, done]
    options= {"cores" : 1, 'memory': "128g", 'walltime': "02:00:00", 'account': "primatediversity"}

    spec = job_header+'''
    mkdir -p "$(dirname {out})"
    python het_stats.py {bcf_file} {ploidy_file} {mask_file} {sexes_file} {out}
    touch {done}
'''.format(bcf_file=bcf_file,ploidy_file=ploidy_file,mask_file=mask_file,
           sexes_file=sexes_file, out=out, done=done)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

