import os,sys
include: 'sample_table.smk'

rule get_all_contigs:
    input:
        expand("{sample}/Viruses/VIBRANT_{sample}_contigs/VIBRANT_phages_{sample}_contigs/{sample}_contigs.phages_combined.fna",
               sample=get_all_samples())
    output:
        temp("viruses/concatenated.fasta")
    shell:
        "cat {input} > {output}"


rule deduplicate:
    input:
        rules.get_all_contigs.output
    output:
        out="viruses/deduplicated.fasta.gz"
    params:
        minoverlap=50,
        maxedits=10,
        findoverlap=True,
        cluster=True,
        processclusters=True,
        command= lambda wc,input: f"cat {input[0]} | dedupe.sh  in=stdin.fasta",
    resources:
        time= 1,
        mem=5
    shadow:
        "minimal"
    log:
        "logs/viruses/deduplicate.log"
    benchmark:
        "logs/benchmark/deduplicate_viruses.txt"
    conda:
        "../envs/bbmap.yaml"
    threads:
        config['threads']
    script:
        "../scripts/runBB.py"



rule sketch:
    input:
        rules.deduplicate.output
    output:
        out="viruses/bbsketch/deduplicated.sketch.gz"
    params:
        k= config['bbsketch']['k'],
        translate=True,
        overwrite=True,
        command= "bbsketch.sh persequence",
    resources:
        time= 10,
        mem= config['mem']
    shadow:
        "minimal"
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
    benchmark:
        "logs/benchmark/bbsketch/alltoall_viruses.txt"
    log:
        "logs/virueses/alltoall.log"
    threads:
        config['threads']
    script:
        "../scripts/runBB.py"


#shHF7/Viruses/VIBRANT_shHF7_contigs/VIBRANT_results_shHF7_contigs/VIBRANT_genome_quality_shHF7_contigs.tsv
# shHF7/Viruses/VIBRANT_shHF7_contigs/VIBRANT_results_shHF7_contigs/VIBRANT_summary_results_shHF7_contigs.tsv
#shHF7/Viruses/VIBRANT_shHF7_contigs/VIBRANT_results_shHF7_contigs/VIBRANT_annotations_shHF7_contigs.tsv
