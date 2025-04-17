rule hisat2_align:
	input:
		reads=["out/trimmed/{sample}_1.fq", "out/trimmed/{sample}_2.fq"],
		idx="out/hisat2_index",
	output:
		"out/mapped/{sample}.bam",
	log:
		"out/logs/hisat2_align/{sample}.log",
	conda:
		"../envs/hisat2.yaml",
	params:
		extra="--dta",
	threads: 8
	script:
		"../scripts/hisat2_align.py"