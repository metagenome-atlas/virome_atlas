import os
import pandas as pd
from snakemake.utils import validate


validate(config, "../config/config.schema.yaml")

DBDIR = config['database_dir']
VIBRANT_DBDIR= os.path.join(DBDIR,'VIBRANT')
VIBRANT_downloaded_flag=os.path.join(VIBRANT_DBDIR,'downloaded_vibrant_data')

include: 'sample_table.smk'


VIBRANT_OUTPUT_CONTIGS= "{sample}/Viruses/VIBRANT_{sample}_contigs/VIBRANT_phages_{sample}_contigs/{sample}_contigs.phages_combined.fna"
VIBRANT_OUTPUT_TABLES= [
    "{sample}/Viruses/VIBRANT_{sample}_contigs/VIBRANT_results_{sample}_contigs/VIBRANT_genome_quality_{sample}_contigs.tsv",
    "{sample}/Viruses/VIBRANT_{sample}_contigs/VIBRANT_results_{sample}_contigs/VIBRANT_summary_results_{sample}_contigs.tsv",
    "{sample}/Viruses/VIBRANT_{sample}_contigs/VIBRANT_results_{sample}_contigs/VIBRANT_annotations_{sample}_contigs.tsv"
]

rule run_vibrant:
    input:
        contigs= "{sample}/{sample}_contigs.fasta",
        database= VIBRANT_downloaded_flag
    output:
        VIBRANT_OUTPUT_CONTIGS,
        VIBRANT_OUTPUT_TABLES
    params:
        folder= "{sample}/Viruses"
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
        plot= "-no_plot" if ~config['vibrant_plot'] else ""
    shell:
        "VIBRANT_run.py -i {input.contigs} -t {threads} -folder {params.folder} "
        " -l {params.min_contig_length} -o {params.minimum_orfs} {params.plot}"



localrules: download_vibrant

rule download_vibrant:
    output:
        touch(VIBRANT_downloaded_flag)
    shadow:
        "minimal"
    conda:
        "../envs/vibrant.yaml"
    threads:
        1
    shell:
        'download-db.sh'
