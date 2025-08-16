rule fastp_pe:
		input:
				get_fq
		output:
				trimmed=["out/trimmed/{sample}_1.fq", "out/trimmed/{sample}_2.fq"],
				html="out/reports/fastp_reports_html/{sample}.html",
				json="out/reports/json/{sample}.json"
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