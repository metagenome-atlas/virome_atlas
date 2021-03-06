

include: "rules/vibrant.smk"
include: "rules/binning.smk"

rule all:
    input:
        expand("{sample}/Viruses",sample=get_all_samples())


rule all_bins:
    input:
        expand("{sample}/binning/viralbins/cluster_attribution.tsv",sample=get_all_samples())


for r in workflow.rules:
    if not "mem" in r.resources:
        r.resources["mem"]=config["mem"]
    if not "time" in r.resources:
        r.resources["time"]=config["time"]


#
