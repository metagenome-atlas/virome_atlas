
def gen_names_for_range(N,prefix='',start=1):
    """generates a range of IDS with leading zeros so sorting will be ok"""
    n_leading_zeros= len(str(N))
    format_int=prefix+'{:0'+str(n_leading_zeros)+'d}'
    return [format_int.format(i) for i in range(start,N+start)]



BINNING_CONTIGS= "{sample}/{sample}_contigs.fasta"


## METABAT
rule get_metabat_depth_file:
    input:
        bam = lambda wc: expand("{sample}/sequence_alignment/{sample}.bam",
                     sample=wc.sample)
    output:
        "{sample}/binning/viralbins/metabat_depth.txt"
    log:
        "{sample}/logs/binning/metabat_depth.txt"
    conda:
        "../envs/metabat.yaml"
    threads:
        config['threads']
    resources:
        mem = config["mem"]
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam} \
            &> {log}
        """


rule metabat:
    input:
        depth_file = rules.get_metabat_depth_file.output,
        contigs = BINNING_CONTIGS
    output:
        "{sample}/binning/viralbins/cluster_attribution.tmp",
    params:
          sensitivity = 2000,
          min_contig_len = 1500,
    benchmark:
        "logs/benchmarks/binning/viralbins/{sample}.txt"
    log:
        "{sample}/logs/binning/viralbins.txt"
    conda:
        "../envs/metabat.yaml"
    threads:
        config["threads"]
    resources:
        mem = config["mem"]
    shell:
        """
        metabat2 -i {input.contigs} \
            --abdFile {input.depth_file} \
            --minContig {params.min_contig_len} \
            --numThreads {threads} \
            --maxEdges {params.sensitivity} \
            --saveCls --noBinOut \
            -o {output} \
            &> {log}
        """


localrules: get_bins
rule get_unique_cluster_attribution:
    input:
        "{sample}/binning/{binner}/cluster_attribution.tmp"
    output:
        "{sample}/binning/{binner}/cluster_attribution.tsv"
    run:
        import pandas as pd
        import numpy as np


        d= pd.read_csv(input[0],index_col=0, squeeze=True, header=None,sep='\t')

        assert type(d) == pd.Series, "expect the input to be a two column file: {}".format(input[0])

        old_cluster_ids = list(d.unique())
        if 0 in old_cluster_ids:
            old_cluster_ids.remove(0)

        map_cluster_ids = dict(zip(old_cluster_ids,
                                   gen_names_for_range(
                                       len(old_cluster_ids),
                                        prefix="{sample}_{binner}_".format(**wildcards)
                                         )
                                   ))

        new_d= d.map(map_cluster_ids)
        new_d.dropna(inplace=True)
        if new_d.shape[0]==0:
            logger.error(f"No bins detected with binner {wildcards.binner} in sample {wildcards.sample}.\n"
                          "This will break the continuationof the pipeline. "
                          "Check what happened. Maybe the the assembly is too small. "
                          "You can either remove the binner (for all samples) from the config.yaml file or the sample from the sample.tsv"
                            )
            raise Exception("No bins detected with binner {wildcards.binner} in sample {wildcards.sample}.")
        new_d.to_csv(output[0],sep='\t')
#


rule get_bins:
    input:
        cluster_attribution = "{sample}/binning/{binner}/cluster_attribution.tsv",
        contigs= BINNING_CONTIGS
    output:
        directory("{sample}/binning/{binner}/bins")
    conda:
        "../envs/sequence_utils.yaml"
    log:
        "{sample}/logs/binning/get_bins_{binner}.log"
    script:
        "../scripts/get_fasta_of_bins.py"
