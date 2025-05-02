"""
------------------------------------------------------------------------------------------------------------------------
This file containes gwf functions of the pipeline.
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Authors: Juraj Bergman, Vasili Pankratov, Bjarke M. Pedersen
Date: 26/03/2025
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
    conda activate /home/bjarkemp/miniforge3/envs/gatksharedenv
    echo "CONDA:" $CONDA_DEFAULT_ENV
'''

default_options = {"cores" : 1, 'memory': "8g", 'walltime': "1-00:00:00", 'account': "primatediversity"}

#### Templates
"""
------------------------------------------------------------------------------------------------------------------------
A. REFERENCE-ASSOCIATED JOBS
------------------------------------------------------------------------------------------------------------------------
"""

def download_ref(ftp, out, done):
    '''Download genbank reference.'''
    inputs = []
    outputs = [done]
    options = default_options.copy()

    spec = job_header+f'''
    rm -f {out}
    curl {ftp} | gzip -d > {out}
    touch {done}
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mask_reference(input, bed, output, done_prev, done):
    """Make reference fasta with contig > 1000 bp in length."""
    inputs = [done_prev]
    outputs = [done]
    options = {"cores" : 1, 'memory': "8g", 'walltime': "1-00:00:00", 'account': "primatediversity"}

    spec = job_header+'''
    bedtools maskfasta -fi {input} -bed {bed} -fo {output}
    touch {done}
    '''.format(input = input, bed = bed, output = output, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def cut_contigs(input, output, min_length, done_prev, done):
    """Make reference fasta with contig > 1000 bp in length."""
    inputs = [done_prev]
    outputs = [done]
    options = {"cores" : 1, 'memory': "8g", 'walltime': "1-00:00:00", 'account': "primatediversity"}

    spec = job_header+'''
    reformat.sh -Xmx8g in={input} out={output} minlength={min_length} overwrite=true
    touch {done}
    '''.format(input = input, output = output, min_length=min_length, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def make_fasta(input, done_prev, done):
    """Make bwa and faidx index."""
    inputs = [done_prev]
    outputs = [done]
    options = {"cores" : 1, 'memory': "16g", 'walltime': "1-00:00:00", 'account': "primatediversity"}

    spec = job_header+'''
    bwa index -a bwtsw {newFasta}

    samtools faidx {newFasta}

    gatk CreateSequenceDictionary -R {newFasta} -O {dict}

    touch {done}
    '''.format(newFasta = input.split(".f")[0] + "_LargerThan1000bp.fasta",
               dict = input.split(".f")[0] + "_LargerThan1000bp.dict", done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def make_regions(refFolder, input, done_prev, done):
    """Make regions file."""
    inputs = [done_prev]
    outputs = [done]
    options = default_options.copy()

    spec = job_header+'''
    python /faststorage/project/primatediversity/people/juraj/python_scripts/make_regions_file.py {refFolder} {chrs}
    touch {done}
    '''.format(refFolder = refFolder, chrs = input.split(".f")[0] + "_LargerThan1000bp.fasta.fai", done=done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

"""
------------------------------------------------------------------------------------------------------------------------
B. SAMPLE DOWNLOAD AND PREPARATION FOR MAPPING
------------------------------------------------------------------------------------------------------------------------
"""

def download_pe2(srr, out, done):
    """Download paired end fastq reads."""
    inputs = []
    outputs = [done]
    options = {'cores': 1, 'memory': '4g', 'walltime': "06:00:00", 'account':"primatediversity"}

    spec = job_header+"""
    rm -f {prev_file}
    wget --progress=dot:giga --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P {out} {srr}
    touch {done}
    """.format(prev_file = out + srr.split("/")[-1], srr = srr, out = out, done = done)

    return AnonymousTarget(inputs = inputs, outputs = outputs, options=options, spec = spec)

def download_per_individual(group, ind, run_accessions, md5s, done):
    """Download paired end fastq reads."""
    inputs = []
    outputs = [done] + ["/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/" + group + "/done/download_pe2_" + group + "_" + run_accession.replace("/", "_") for run_accession in run_accessions.split(",")]
    options = {'cores': 1, 'memory': '8g', 'walltime': "06:00:00", 'account':"primatediversity"}

    spec = job_header+f"""
    python /faststorage/project/primatediversity/people/juraj/python_scripts/download_fastq_per_individual_check_md5.py {group} {ind} {run_accessions} {md5s}
    """

    return AnonymousTarget(inputs = inputs, outputs = outputs, options=options, spec = spec)

def prefetch(srr, out, done):
    """Download and split paired end fastq reads."""
    inputs = []
    outputs = [done]
    options = {'cores': 1, 'memory': '4g', 'walltime': "06:00:00", 'account':"primatediversity"}

    spec = job_header+"""
    cd {out}
    prefetch --max-size 1t   {srr}
    touch {done}
    """.format(srr = srr, out = out, done = done)

    return AnonymousTarget(inputs = inputs, outputs = outputs, options=options, spec = spec)

def split_srr(srr, out, done_prev, done):
    """Split paired end fastq reads."""
    inputs = [done_prev]
    outputs = [done+"_1.fastq.gz", done+"_2.fastq.gz"]
    options = {'cores': 1, 'memory': '4g', 'walltime': "12:00:00", 'account':"primatediversity"}

    spec = job_header+"""
    conda activate workflowMap
    fastq-dump --gzip --split-files -O {out} {srr}
    rm {srr}
    touch {done1}
    touch {done2}
    """.format(srr = out+srr, out = out, done1 = done+"_1.fastq.gz", done2 = done+"_2.fastq.gz")

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def concatfastqs(group, ind, fastqs_1, fastqs_2, prev_done, done):
    '''
    Concatenate fastqs
    '''
    inputs    = prev_done
    outputs   = [done]
    options   = default_options.copy()
    spec      = job_header+f'''
    out_dir=/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq

    mkdir -p ${{out_dir}}

    cat {" ".join(fastqs_1)} > ${{out_dir}}/{ind}_R1.fastq.gz
    cat {" ".join(fastqs_2)} > ${{out_dir}}/{ind}_R2.fastq.gz
    
    rm {" ".join(fastqs_1)}
    rm {" ".join(fastqs_2)}
    
    touch {done}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def renamefastqs(group, ind, fastq_1, fastq_2, prev_done, done):
    '''
    Rename fastqs
    '''
    inputs    = prev_done
    outputs   = [done]
    options   = default_options.copy()
    spec      = job_header+f'''
    out_dir=/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/fastq

    mkdir -p ${{out_dir}}

    mv {fastq_1} ${{out_dir}}/{ind}_R1.fastq.gz
    mv {fastq_2} ${{out_dir}}/{ind}_R2.fastq.gz

    touch {done}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def makeuBAM(group, ind, fastq1, fastq2, prev_done, done):
    '''
    Make uBAM file
    '''
    inputs  = prev_done
    outputs = [done]
    options = {"cores" : 1, 'memory': "8g", 'walltime': "1-00:00:00", 'account': "primatediversity"}
    spec    = job_header+f"""
    out_dir=/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam

    mkdir -p ${{out_dir}}

    picard FastqToSam --FASTQ  {fastq1}                                          \
                      --FASTQ2 {fastq2}                                          \
                      --OUTPUT /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/u{ind}.bam \
                      --SAMPLE_NAME {ind}                                        \
                      #--IO_SIZE null                                             \
                      --TMP_DIR /scratch/$SLURM_JOB_ID
    rm {fastq1} {fastq2}   
    touch {done}
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def splituBAM(group, ind, prev_done, done):
    '''
    Split uBAMs
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec    = job_header+f"""
    splitubamdir=/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}
    rm    -fr ${{splitubamdir}}
    mkdir -p  ${{splitubamdir}} 

    picard SplitSamByNumberOfReads -I /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/u{ind}.bam \
                                   -O ${{splitubamdir}}                \
                                   --CREATE_INDEX true                 \
                                   -N_READS 48000000

    ls ${{splitubamdir}} | wc -l > /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_nsplitubams.txt
    
    rm /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/u{ind}.bam
    
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

"""
------------------------------------------------------------------------------------------------------------------------
C. MAPPING
------------------------------------------------------------------------------------------------------------------------
"""

def shardstr(ishard):
    return (4-len(str(ishard+1)))*"0"+str(ishard+1)

def markadapt(group, ind, shard, prev_done, done):
    '''
    Mark adaptor sequences in uBAM file
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    picard MarkIlluminaAdapters -I /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}/shard_{shard}.bam           \
                                -O /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}/shard_{shard}_markadapt.bam \
                                -M /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}/shard_{shard}_markadapt.txt \
                                --TMP_DIR /scratch/$SLURM_JOB_ID

    rm -f /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/u{ind}.bam
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mapBAM(group, ind, shard, ref, prev_done, done):
    '''
    Mapping bam file
    '''
    inputs  = prev_done
    outputs = [done]
    # options = default_options.copy()
    options = {"cores" : 1, 'memory': "16g", 'walltime': "1-00:00:00", 'account': "primatediversity"}
    spec = job_header+f"""
    picard SamToFastq -I           /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}/shard_{shard}_markadapt.bam       \
                      --FASTQ      /dev/stdout                                                        \
                      --INTERLEAVE true                                                               \
                      --TMP_DIR    /scratch/$SLURM_JOB_ID                                             \
        | bwa mem -M                    \
                  -t {options["cores"]} \
                  -p {ref}              \
                  /dev/stdin            \
        | picard MergeBamAlignment --ALIGNED_BAM                  /dev/stdin                                                          \
                                   --UNMAPPED_BAM                 /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}/shard_{shard}_markadapt.bam        \
                                   --OUTPUT                       /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}/shard_{shard}_markadapt_mapped.bam \
                                   -R                             {ref}                                                               \
                                   --CREATE_INDEX                 true                                                                \
                                   --INCLUDE_SECONDARY_ALIGNMENTS false                                                               \
                                   --TMP_DIR                      /scratch/$SLURM_JOB_ID

    rm -f /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}/shard_{shard}.bam
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mergeBAMs(group, ind, bams, prev_done, done):
    '''
    Merge bam files.
    '''
    inputs  = prev_done
    outputs = [done]
    # options = default_options.copy()
    options = {"cores": 1, 'memory': "16g", 'walltime': "1-00:00:00", 'account': "primatediversity"}
    spec = job_header+f"""
    picard MergeSamFiles {" ".join([f"-I {bam} " for bam in bams])}                  \
                         -O /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged.bam  \
                         --SORT_ORDER   queryname                                    \
                         --TMP_DIR      /scratch/$SLURM_JOB_ID

    rm -f /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}/shard_*_markadapt.bam
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def markduplicates(group, ind, prev_done, done):
    '''
    Mark and remove duplicates
    '''
    inputs  = prev_done
    outputs = [done]
    # options = default_options.copy()
    options = {"cores": 1, 'memory': "16g", 'walltime': "12:00:00", 'account': "primatediversity"}
    spec = job_header+f"""
    picard MarkDuplicates -I /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged.bam                \
                          -M /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged_markduplicates.txt \
                          -O /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged_markduplicates.bam \
                          --REMOVE_DUPLICATES true                                                  \
                          --CREATE_INDEX true                                                       \
                          --TMP_DIR /scratch/$SLURM_JOB_ID

    rm -fr /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/split_uBAM{ind}
    rm -f /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged_markduplicates.txt
    rm /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged.bam
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def coordsort(group, ind, prev_done, done):
    '''
    Sort BAM by coordinates
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    picard SortSam -I /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged_markduplicates.bam           \
                   -O /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bam \
                   -SO coordinate                                                                      \
                   --CREATE_INDEX true                                                                 \
                   --TMP_DIR /scratch/$SLURM_JOB_ID

    rm -f /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged_markduplicates.bam
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def cov(group, ind, regions, chromosomes, starts, ends, prev_done, done):
    '''
    Get average coverage across chromosomes at covered sites, for files already coordinate-sorted
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    covdir=/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/cov
    # rm    -fr ${{covdir}}
    mkdir -p  ${{covdir}} 

    regions=({" ".join(regions)})
    chromosomes=({" ".join(chromosomes)})
    starts=({" ".join([str(s) for s in starts])})
    ends=({" ".join([str(e) for e in ends])})
    length=${{#regions[@]}}

    # Iterate over both arrays simultaneously
    for ((i = 0; i < length; i++)); do
        region=${{regions[i]}}
        chrom=${{chromosomes[i]}}
        start=${{starts[i]}}
        end=${{ends[i]}}

        echo ${{region}} ${{chrom}}
        date

        samtools depth -r ${{chrom}}:${{start}}-${{end}} /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/bam/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bam \
            | awk '{{sum += $3}} END {{if(sum == 0 || NR == 0){{cov=0}}else{{cov=sum/NR}};print "'${{region}}'\t'${{chrom}}'\t'${{start}}'\t'${{end}}'\t"NR"\t"sum"\t"cov}}' >> ${{covdir}}/{ind}.cov

    done

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

"""
------------------------------------------------------------------------------------------------------------------------
D. CALLING AND GENOTYPING
------------------------------------------------------------------------------------------------------------------------
"""
def find_chrX(subset_file, ref_folder, contigs, minlen, prev_done, done):
    inputs = prev_done
    outputs = [done]
    options = default_options.copy()

    spec = job_header + f'''
    python /faststorage/project/primatediversity/people/juraj/python_scripts/find_chrX.py -s {subset_file} -r {ref_folder} -c {contigs} -m {minlen}
  
    touch {done}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def make_simplified_batch_file(regions_file, ref_folder, prev_done, done):
    inputs = prev_done
    outputs = [done]
    options = default_options.copy()

    spec = job_header + f'''
    python /faststorage/project/primatediversity/people/juraj/python_scripts/make_simplified_batch_file.py {regions_file} {ref_folder}
  
    touch {done}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def call_batch(group, ind, batch, chromosomes, starts, ends, ref, ploidy, prev_done, done):
    '''
    Call a batch of chromosomes per individual
    '''
    inputs = prev_done
    outputs = [done]
    # options = default_options.copy()
    options = {"cores": 1, 'memory': "8g", 'walltime': "12:00:00", 'account': "primatediversity"}

    spec = job_header + """
            echo -e "{bed}" > /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}.bed

            samtools view -b -L /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}.bed {dir}bam/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bam  >  /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}.bam

            gatk AddOrReplaceReadGroups -I             /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}.bam      \
                                        -O             /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}_readgroup.bam \
                                        -LB            lib1                                                   \
                                        -PL            ILLUMINA                                               \
                                        -PU            unit1                                                  \
                                        -SM            {ind}                                                  \
                                        --CREATE_INDEX true                                                   \
                                        --TMP_DIR      /scratch/$SLURM_JOB_ID

            gatk HaplotypeCaller -R                       {ref}                                                  \
                                 -I                       /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}_readgroup.bam \
                                 -L {regions}  \
                                 -ploidy                   {ploidy}                                               \
                                 --native-pair-hmm-threads {cores}                                                \
                                 -ERC                      BP_RESOLUTION                                          \
                                 -O                        {dir}gVCF/{ind}_batch_{batch}_ploidy_{ploidy}.gvcf.gz 

            gatk IndexFeatureFile -I {dir}/gVCF/{ind}_batch_{batch}_ploidy_{ploidy}.gvcf.gz

            touch {done}
            """.format(bed="\n".join(chromosomes[j] + "\t" + str(starts[j] - 1) + "\t" + str(ends[j]) for j in range(len(chromosomes))),
                       group=group, ind=ind, batch=batch, ref=ref, ploidy=ploidy, cores=options["cores"],
                       regions=" -L ".join(chromosomes[j] + ":" + str(starts[j]) + "-" + str(ends[j]) for j in range(len(chromosomes))),
                       dir="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/" + group + "/",
                       done=done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_batch_metadata(loc, group, regions_file, inds, sexes, ref_folder, prev_done, done):
    inputs = prev_done
    outputs = [done]
    options = {"cores": 1, 'memory': "8g", 'walltime': "00:05:00", 'account': "primatediversity"}

    spec = job_header + f'''
    python /faststorage/project/primatediversity/people/juraj/python_scripts/make_loc_metadata.py {loc} {group} {regions_file} {inds} {sexes} {ref_folder}
  
    touch {done}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def call_batch_with_bed(group, ind, batch, bed, intervals, ref, ploidy, prev_done, done):
    '''
    Call a batch of chromosomes per individual
    '''
    inputs = prev_done
    outputs = [done]
    # options = default_options.copy()
    options = {"cores": 1, 'memory': "8g", 'walltime': "24:00:00", 'account': "primatediversity"}

    spec = job_header + """

            samtools view -b -L {bed} {dir}bam/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bam  >  /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}.bam

            gatk AddOrReplaceReadGroups -I             /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}.bam      \
                                        -O             /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}_readgroup.bam \
                                        -LB            lib1                                                   \
                                        -PL            ILLUMINA                                               \
                                        -PU            unit1                                                  \
                                        -SM            {ind}                                                  \
                                        --CREATE_INDEX true                                                   \
                                        --TMP_DIR      /scratch/$SLURM_JOB_ID

            gatk HaplotypeCaller -R                       {ref}                                                  \
                                 -I                       /scratch/$SLURM_JOB_ID/{ind}_batch_{batch}_ploidy_{ploidy}_readgroup.bam \
                                 -L {intervals}  \
                                 -ploidy                   {ploidy}                                               \
                                 --native-pair-hmm-threads {cores}                                                \
                                 -ERC                      BP_RESOLUTION                                          \
                                 -O                        {dir}gVCF/{ind}_batch_{batch}_ploidy_{ploidy}.gvcf.gz 

            gatk IndexFeatureFile -I {dir}/gVCF/{ind}_batch_{batch}_ploidy_{ploidy}.gvcf.gz

            touch {done}
            """.format(bed=bed,
                       group=group, ind=ind, batch=batch, ref=ref, ploidy=ploidy, cores=options["cores"],
                       intervals=intervals,
                       dir="/faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/" + group + "/",
                       done=done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def make_genDB_folder_and_map(group, batch, inds, ploidies, folder_name, prev_done, done):
    inputs = prev_done
    outputs = [done]
    options = {"cores": 1, 'memory': "8g", 'walltime': "00:05:00", 'account': "primatediversity"}

    spec = job_header + '''
    python /faststorage/project/primatediversity/people/juraj/python_scripts/make_genDB_folder_and_map.py {group} {batch} {inds} {ploidies} {folder_name}
  
    touch {done}
    '''.format(group=group, batch=batch, inds=inds, ploidies=ploidies, folder_name=folder_name, done=done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_genDB(group, batch, fploidy, mploidy, chromosomes, starts, ends, prev_done, done):
    '''
    Create a GenomicsDB
    '''
    inputs = [prev_done]
    outputs = [done]
    options = default_options.copy()
    spec = job_header + """
    gatk GenomicsDBImport --sample-name-map                         /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/cohort.sample_map \
                          --genomicsdb-workspace-path               /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}   \
                          --tmp-dir                                 /scratch/$SLURM_JOB_ID      \
                          --genomicsdb-shared-posixfs-optimizations true                        \
                          --genomicsdb-vcf-buffer-size              4194304                     \
                          --reader-threads                          {cores} \
                          -L {intervals}


    chmod -R 777 /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    rm -rf /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    mv /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    touch {done}
    """.format(group=group, batch=batch, fploidy = fploidy, mploidy = mploidy, done=done, cores=options["cores"],
               intervals=" -L ".join(chromosomes[j] + ":" + str(starts[j]) + "-" + str(ends[j]) for j in range(len(chromosomes))))
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def make_genDB_with_bed(group, batch, fploidy, mploidy, intervals, prev_done, done):
    '''
    Create a GenomicsDB
    '''
    inputs = [prev_done]
    outputs = [done]
    options = default_options.copy()
    spec = job_header + """
    gatk GenomicsDBImport --sample-name-map                         /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/cohort.sample_map \
                          --genomicsdb-workspace-path               /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}   \
                          --tmp-dir                                 /scratch/$SLURM_JOB_ID      \
                          --genomicsdb-shared-posixfs-optimizations true                        \
                          --genomicsdb-vcf-buffer-size              4194304                     \
                          --reader-threads                          {cores} \
                          -L {intervals}


    chmod -R 777 /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    rm -rf /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    mv /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    touch {done}
    """.format(group=group, batch=batch, fploidy = fploidy, mploidy = mploidy, done=done, cores=options["cores"],
               intervals=intervals)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def GenotypeGVCFs_subbatch(group, batch, fploidy, mploidy, subbatch, chromosome, start, end, out, ref, cores, memory, prev_done, done):
    '''
    Genotype.
    '''
    inputs = prev_done
    outputs = [done]
    options = {"cores": cores, 'memory': memory, 'walltime': "1-00:00:00", 'account': "primatediversity"}
    spec = job_header + """
    cp -r /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    gatk GenotypeGVCFs --include-non-variant-sites                                                                \
                       -R                                 {ref}                                                   \
                       --variant                          gendb:///scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} \
                       -O                                 /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_subbatch_{subbatch}_gt.gvcf.gz \
                       --max-genotype-count               3000                                                    \
                       --genomicsdb-max-alternate-alleles 3001                                                    \
                       --max-alternate-alleles            6 \
                       -L {interval} 

    chmod -R 777 /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_subbatch_{subbatch}_gt.gvcf.gz
    mv /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_subbatch_{subbatch}_gt.gvcf.gz {out}

    touch {done}
    """.format(group=group, batch=batch, fploidy=fploidy, mploidy=mploidy, subbatch=subbatch,  ref=ref, done = done,  out = out,
               interval = chromosome + ":"  + str(start) + "-" + str(end))
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def GenotypeGVCFs_subbatch_with_bed(group, batch, fploidy, mploidy, subbatch, interval, out, ref, cores, memory, prev_done, done):
    '''
    Genotype.
    '''
    inputs = prev_done
    outputs = [done]
    options = {"cores": cores, 'memory': memory, 'walltime': "1-00:00:00", 'account': "primatediversity"}
    spec = job_header + """
    cp -r /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    gatk GenotypeGVCFs --include-non-variant-sites                                                                \
                       -R                                 {ref}                                                   \
                       --variant                          gendb:///scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} \
                       -O                                 /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_subbatch_{subbatch}_gt.gvcf.gz \
                       --max-genotype-count               3000                                                    \
                       --genomicsdb-max-alternate-alleles 3001                                                    \
                       --max-alternate-alleles            6 \
                       -L {interval} 

    chmod -R 777 /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_subbatch_{subbatch}_gt.gvcf.gz
    mv /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_subbatch_{subbatch}_gt.gvcf.gz {out}

    touch {done}
    """.format(group=group, batch=batch, fploidy=fploidy, mploidy=mploidy, subbatch=subbatch,  ref=ref, done = done,  out = out,
               interval = interval)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def GenotypeGVCFs(group, batch, fploidy, mploidy, chromosomes, starts, ends, out, ref, cores, memory, prev_done, done):
    '''
    Genotype
    '''
    inputs = prev_done
    outputs = [done]
    options = {"cores": cores, 'memory': memory, 'walltime': "1-00:00:00", 'account': "primatediversity"}
    spec = job_header + """
    cp -r /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    gatk GenotypeGVCFs --include-non-variant-sites                                                                \
                       -R                                 {ref}                                                   \
                       --variant                          gendb:///scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} \
                       -O                                 /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_gt.gvcf.gz \
                       -L {intervals} 

    chmod -R 777 /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_gt.gvcf.gz
    mv /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_gt.gvcf.gz {out}

    touch {done}
    """.format(group=group, batch=batch, fploidy=fploidy, mploidy=mploidy, ref=ref, done=done, out = out,
               intervals=" -L ".join(chromosomes[j] + ":" + str(starts[j]) + "-" + str(ends[j]) for j in range(len(chromosomes))))
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def GenotypeGVCFs_with_bed(group, batch, fploidy, mploidy, intervals, out, ref, cores, memory, prev_done, done):
    '''
    Genotype
    '''
    inputs = prev_done
    outputs = [done]
    options = {"cores": cores, 'memory': memory, 'walltime': "1-00:00:00", 'account': "primatediversity"}
    spec = job_header + """
    cp -r /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} /scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    gatk GenotypeGVCFs --include-non-variant-sites                                                                \
                       -R                                 {ref}                                                   \
                       --variant                          gendb:///scratch/$SLURM_JOB_ID/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy} \
                       -O                                 /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_gt.gvcf.gz \
                       -L {intervals} 

    chmod -R 777 /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_gt.gvcf.gz
    mv /scratch/$SLURM_JOB_ID/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_gt.gvcf.gz {out}

    touch {done}
    """.format(group=group, batch=batch, fploidy=fploidy, mploidy=mploidy, ref=ref, done=done, out = out,
               intervals=intervals)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def bcftoolsconcat(vcfs, vcf, prev_done, done):
    '''
    concatenate vcfs
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    bcftools concat {vcfs} -O z -o {vcf}

    rm {vcfs}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def picardconcat(vcfs, vcf, prev_done, done):
    '''
    concatenate  and sort vcfs
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+"""

    picard SortVcf -I {vcfs} -O {vcf}

    rm -f {rmvcfs}

    touch {done}
    """.format(vcfs = " -I ".join(vcfs), vcf = vcf, rmvcfs = " ".join(vcfs), done = done)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def renameGVCF(vcf_in, vcf_out, prev_done, done):
    '''
    Rename GVCF.
    '''
    inputs    = prev_done
    outputs   = [done]
    options   = default_options.copy()
    spec      = job_header+f'''

    mv {vcf_in} {vcf_out}

    touch {done}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def IndexGVCFs(group, batch, fploidy, mploidy, prev_done, done):
    '''
    Index GVCF.
    '''

    inputs = [prev_done]
    outputs = [done]
    options = default_options.copy()

    spec = job_header + f"""
    gatk IndexFeatureFile -I /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/gVCF/{group}_batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}_gt.gvcf.gz

    rm -r /faststorage/project/primatediversity/data/gVCFs_recalling_10_12_2024/{group}/GenomicsDB/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}/batch_{batch}_fploidy_{fploidy}_mploidy_{mploidy}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

"""
------------------------------------------------------------------------------------------------------------------------
E. CLEAN UP
------------------------------------------------------------------------------------------------------------------------
"""

def bam_to_cram(path_to_bam_folder, ind, ref, done):
    """Convert bam to cram."""
    inputs = []
    outputs = [done]
    options = {'cores': 1, 'memory': "8g", 'walltime': "12:00:00", 'account':"primatediversity"}

    spec = f"""
    samtools view -C -T {ref} -o {path_to_bam_folder}/{ind}_markadapt_mapped_merged_markduplicates_coordsort.cram {path_to_bam_folder}/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bam
    samtools index {path_to_bam_folder}/{ind}_markadapt_mapped_merged_markduplicates_coordsort.cram
    
    rm -f {path_to_bam_folder}/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bam
    rm -f {path_to_bam_folder}/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bai

    rm -f {path_to_bam_folder}/u{ind}.bam
    rm -f {path_to_bam_folder}/{ind}_nsplitubams.txt

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def removeIndGVCFs(cohort_map, prev_done, done):
    '''
    Remove un-genotyped individual-based GVCFs.
    '''
        
    inputs = [prev_done]
    outputs = [done]
    options = {"cores": 1, 'memory': "8g", 'walltime': "00:05:00", 'account': "primatediversity"}

    spec = job_header + f'''
    python /faststorage/project/primatediversity/people/juraj/python_scripts/removeIndGVCFs.py {cohort_map}
    touch {done}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)