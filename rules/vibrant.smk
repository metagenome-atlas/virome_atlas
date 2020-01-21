import os
import pandas as pd
from snakemake.utils import validate


validate(config, "../config/config.schema.yaml")
print(config)

DBDIR = config['database_dir']
VIBRANT_DBDIR= os.path.join(DBDIR,'VIBRANT')


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

    SampleTable = pd.read_csv(SampleTableFile, sep='\t').set_index("Samples", drop=False)

    validate(SampleTable, "../config/samples.schema.yaml")

    return list(SampleTable.index)



rule vibrant:
    input:
        expand("{sample}/Viruses",sample=get_all_samples())


rule run_vibrant:
    input:
        contigs= "{sample}/{sample}_contigs.fasta",
        database= VIBRANT_DBDIR
    output:
        directory("{sample}/Viruses")
    resources:
        mem= config.get('mem',70),
        time = config.get('time',24)
    threads:
        config.get("threads",8)
    conda:
        "../vibrant.yaml"
    params:
        min_contig_length=config['vibrant_min_contig_length'],
        minimum_orfs= config['vibrant_minimum_orfs'],
        plot= "-no_plot" if ~config['vibrant_plot'] else ""
    shell:
        'export VIBRANT_DATA_PATH="{input.database}" ; '
        "VIBRANT_run.py -i {input.contigs} -t {threads} -folder {output} "
        " -l {params.min_contig_length} -o {params.minimum_orfs} {params.plot}"



localrules: download_vibrant

rule download_vibrant:
    output:
        VIBRANT_DBDIR
    shadow:
        "minimal"
    conda:
        "../vibrant.yaml"
    threads:
        1
    shell:
        'export VIBRANT_DATA_PATH="{output}"  ; download-db.sh'
