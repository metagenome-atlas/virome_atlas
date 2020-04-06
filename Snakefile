import os

configfile: os.path.join(os.path.dirname(workflow.snakefile),"config/default_config.yaml")
include: "rules/vibrant.smk"
include: "rules/dereplicate.smk"

rule all:
    input:
        expand("{sample}/Viruses",sample=get_all_samples())


rule all_bins:
    input:
        expand(VIBRANT_OUTPUT_CONTIGS,sample=get_all_samples()),
        expand(VIBRANT_OUTPUT_TABLES,sample=get_all_samples())

rule dereplicate:
    input:
        "viruses/bbsketch/dists.tsv"


for r in workflow.rules:
    if not "mem" in r.resources:
        r.resources["mem"]=config["mem"]
    if not "time" in r.resources:
        r.resources["time"]=config["time"]


#
