rule counts:
    input:
        bam="out/mapped/{sample}.bam",
        gtf=config["ref"]["annotation"]
    params:    
        attribute_type=config["featureCounts"]["attribute_type"],
        feature_type=config["featureCounts"]["feature_type"]
    output:
        "out/featurecounts/{sample}.txt"
    log:
        "out/logs/featurecounts/{sample}.log",
    threads:5
    conda:
        "../envs/featurecounts.yaml"
    shell:
        "featureCounts -p \
        -t {params.feature_type} \
        -g {params.attribute_type} \
        -a {input.gtf} \
        -o {output} \
        -T {threads} \
        {input.bam} \
        1>{log} 2>&1"

rule counts_merge:
    input:
        expand("out/featurecounts/{sample}.txt",sample=SAMPLES)
    output:
        merged="out/counts/samples_merged_counts.csv",
        tpm="out/counts/samples_merged_tpm.csv",
    params:
        input_dir="out/featurecounts"
    conda:
        "../envs/featurecounts.yaml"
    script:
        "../scripts/merge_matrix.py"