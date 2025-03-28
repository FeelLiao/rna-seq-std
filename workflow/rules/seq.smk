rule fastp_pe:
    input:
        get_fq
    output:
        trimmed=["out/trimmed/{sample}_1.fq", "out/trimmed/{sample}_2.fq"],
        html="out/trimmed_reports/html/{sample}.html",
        json="out/trimmed_reports/json/{sample}.json"
    log:
        "out/logs/fastp/{sample}.log"
    conda:
        "../envs/fastp.yaml",
    params:
        adapters="--detect_adapter_for_pe",
        extra="",
    threads: 5
    script:
        "../scripts/fastp.py"

rule samtools_flagstat:
    input:
        "out/mapped/{sample}.bam",
    output:
        "out/flagstat/{sample}.txt",
    params:
        extra="",  # optional params string
    threads: 5
    conda:
        "../envs/hisat2.yaml",
    script:
        "../scripts/samtools_flagstat.py"