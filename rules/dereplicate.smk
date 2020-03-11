import os,sys
include: 'sample_table.smk'




rule sketch:
    input:
        expand("{sample}/viruses/.../{sample}_contigs.phages_combined.fna",sample=get_all_samples())
    output:
        out="viruses/bbsketch/combined.sketch.gz"
    params:
        k= config['bbsketch']['k'],
        translate=True,
        overwrite=True,
        command= "bbsketch.sh persequence",
    resources:
        time= 5
    log:
        "logs/viruses/sketch.log"
    benchmark:
        "logs/benchmark/bbsketch/viruses.txt"
    conda:
        "../envs/bbmap.yaml"
    threads:
        config['threads']
    script:
        "../scripts/runBB.py"


rule allvall:
    input:
        ref=rules.sketch.output
    output:
        out="viruses/bbsketch/dists.tsv"
    params:
        amino=True,
        overwrite=True,
        command="comparesketch.sh alltoall",
        prealloc=0.75,
        format=3,
        k=config['bbsketch']['k'],
    shadow:
        "minimal"
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem= config['mem']['large']
    benchmark:
        "logs/benchmark/bbsketch/alltoall_viruses.txt"
    log:
        "logs/virueses/alltoall.log"
    threads:
        config['threads']
    script:
        "../scripts/runBB.py"
