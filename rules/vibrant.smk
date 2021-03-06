import os
import pandas as pd
from snakemake.utils import validate


validate(config, "../config/config.schema.yaml")

DBDIR = config['database_dir']
VIBRANT_DBDIR= os.path.join(DBDIR,'VIBRANT')
VIBRANT_downloaded_flag=os.path.join(VIBRANT_DBDIR,'downloaded_vibrant_data')

def get_all_samples():

    SampleTableFile= config["sampletable"]

    if not os.path.exists(SampleTableFile):
        logger.error(f"The configuration says I have to look for SampleTable at {SampleTableFile}\n"
                      "But this file doesn't exists! The SampleTable should be a table in the format: \n"
                      "Samples\tOPTIONAL_HEADER\n"
                      "Sample1\t\n"
                      "Sample2\t\n"
                      )
        raise IOError("sampletable doesn't exist")

    SampleTable = pd.read_csv(SampleTableFile, sep='\t', index_col=0)

    return list(SampleTable.index)




rule run_vibrant:
    input:
        contigs= "{sample}/{sample}_contigs.fasta",
        database= VIBRANT_downloaded_flag
    output:
        directory("{sample}/Viruses")
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
        "VIBRANT_run.py -i {input.contigs} -t {threads} -folder {output} "
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
