rule post_process:
    input:
        get_final_output(),
    output:
        "out/trimmed_reports/trim_report.csv",
        "out/hisat2_align_report.csv"
    conda:
        "../envs/reports.yaml",
    params:
        trimdrp ="out/trimmed_reports",  
        hisat2log = "out/logs/hisat2_align",
        clean = config["clean"],
    script:
        "../report/post_process.py"

rule reports:
    input:
        "out/trimmed_reports/trim_report.csv",
        "out/hisat2_align_report.csv"
    output:
        "out/reports/report.html",
    conda:
        "../envs/reports.yaml",
    script:
        "../report/reports.Rmd"