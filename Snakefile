

include: "rules/vibrant.smk"
include: "rules/binning.smk"

rule all:
    input:
        expand("{sample}/Viruses",sample=get_all_samples())


rule all_bins:
    input:
        expand("{sample}/binning/viralbins/cluster_attribution.tsv",sample=get_all_samples())
