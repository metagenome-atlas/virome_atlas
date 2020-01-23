

include: "rules/vibrant.smk"

rule all:
    input:
        expand("{sample}/VIBRANT_{sample}_contigs",sample=get_all_samples())
