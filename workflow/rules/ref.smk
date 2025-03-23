rule hisat2_index:
    input:
        fasta=config["ref"]["genome"],
    output:
        directory("out/hisat2_index"),
    params:
        extra="--seed 42",
    conda:
        "../envs/hisat2.yaml",
    log:
        "out/logs/hisat2_index.log",
    threads: 5
    script:
        "../scripts/hisat2_index.py"

#TODO: redirect log to logs/{accession}.log
# try to define log path in params. 
rule download_sra:
    input:
        config["SRA"]["acc_list"]
    output:
        directory(config["SRA"]["output_dir"]),
    conda:
        "../envs/sratools.yaml"
    log:
        "out/logs/sra_download.log",
    params:
        extra="--skip-technical",
    threads: 8  # defaults to 6
    script:
        "../scripts/sra_download.py"