

include: "rules/vibrant.smk"

rule all:
    input:
        expand("{sample}/Viruses",sample=get_all_samples())
