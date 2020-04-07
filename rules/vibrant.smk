import os
import pandas as pd
from snakemake.utils import validate


validate(config, "../config/config.schema.yaml")

DBDIR = config['database_dir']
VIBRANT_DBDIR= os.path.join(DBDIR,'VIBRANT')
VIBRANT_downloaded_flag=os.path.join(VIBRANT_DBDIR,'downloaded_vibrant_data')

include: 'sample_table.smk'


VIBRANT_OUTPUT_CONTIGS= "{sample}/Viruses/VIBRANT_phages_{sample}_contigs/{sample}_contigs.phages_combined.fna"
VIBRANT_OUTPUT_TABLES= [
    "{sample}/Viruses/VIBRANT_results_{sample}_contigs/VIBRANT_genome_quality_{sample}_contigs.tsv",
    "{sample}/Viruses/VIBRANT_results_{sample}_contigs/VIBRANT_summary_results_{sample}_contigs.tsv",
    "{sample}/Viruses/VIBRANT_results_{sample}_contigs/VIBRANT_annotations_{sample}_contigs.tsv"
]

rule run_vibrant:
    input:
        contigs= "{sample}/{sample}_contigs.fasta",
        database= ancient(VIBRANT_downloaded_flag)
    output:
        VIBRANT_OUTPUT_CONTIGS,
        VIBRANT_OUTPUT_TABLES
    resources:
        mem= config.get('mem',70),
        time = config.get('time',24)
    threads:
        config.get("threads",8)
    conda:
        "../envs/vibrant.yaml"
    params:
        min_contig_length=config['vibrant_min_contig_length'],
        minimum_orfs= config['vibrant_minimum_orfs'],
        plot= "-no_plot" if ~config['vibrant_plot'] else "",
        folder= "{sample}/Viruses"
    shell:
        "VIBRANT_run.py -i {input.contigs} -t {threads} -folder {params.folder} "
        " -l {params.min_contig_length} -o {params.minimum_orfs} {params.plot}"

# takeover from old version
localrules: symlink
ruleorder: symlink > run_vibrant
rule symlink:
    input:
        "{sample}/Viruses/VIBRANT_{sample}_contigs/VIBRANT_phages_{sample}_contigs/{sample}_contigs.phages_combined.fna",
        "{sample}/Viruses/VIBRANT_{sample}_contigs/VIBRANT_results_{sample}_contigs/VIBRANT_genome_quality_{sample}_contigs.tsv",
        "{sample}/Viruses/VIBRANT_{sample}_contigs/VIBRANT_results_{sample}_contigs/VIBRANT_summary_results_{sample}_contigs.tsv",
        "{sample}/Viruses/VIBRANT_{sample}_contigs/VIBRANT_results_{sample}_contigs/VIBRANT_annotations_{sample}_contigs.tsv"
    output:
        VIBRANT_OUTPUT_CONTIGS,
        VIBRANT_OUTPUT_TABLES
    run:
        from common import io

        io.symlink_relative([f"VIBRANT_phages_{wildcards.sample}_contigs",
                             f"VIBRANT_results_{wildcards.sample}_contigs"
                             ],
                            input_dir = f"{wildcards.sample}/Viruses/VIBRANT_{wildcards.sample}_contigs",
                            output_dir = f"{wildcards.sample}/Viruses"
                            )



localrules: download_vibrant

rule download_vibrant:
    output:
        touch(VIBRANT_downloaded_flag)
    conda:
        "../envs/vibrant.yaml"
    threads:
        1
    shell:
        'download-db.sh'
