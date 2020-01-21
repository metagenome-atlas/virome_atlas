
import os

def get_and_verify_config_attribute(config_attribute):

    if not config_attribute in config:
        raise IOError(f"'{config_attribute}' is expected to be in the config.yaml file or"
                       "suplied by the '--config' comand line parameter ")

    return config[config_attribute]

DBDIR = get_and_verify_config_attribute('database_dir'])


def get_samples():

    SampleTableFile = get_and_verify_config_attribute('sampletable')

    if not os.path.exists(SampleTableFile):
        raise IOError(f"The configuration says I have to look for SampleTable at {SampleTableFile}\n"
                      "But this file doesn't exists! The SampleTable should be a table in the format: \n"
                      "SAMPLES\tOPTIONAL_HEADER\n"
                      "Sample1\t\n"
                      "Sample2\t\n"
                      )


    SAMPLES = pd.read_csv(SampleTableFile, index_col=0, sep='\t').index

    return list(SAMPLES)



VIBRANT_DBDIR= os.path.join(DBDIR,'VIBRANT')




rule vibrant:
    input:
        expand("{sample}/Viruses",sample=get_samples())

run_vibrant:
    input:
        contigs= "{sample}/{sample}_contigs.fasta"
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
        min_contig_length=3000,
        minimum_orfs= 4,
        plot= "-no_plot" if False else ""
    shell:
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
        "mkdir {output} ; cd {output} ; download-db.sh"
