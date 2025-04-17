# Initial assembly of transcripts with StringTie for each sample 
rule stringtie_initial:
	input: 
		sbam="out/mapped/{sample}.bam",
		anno=config["ref"]["annotation"],
	output:
		"out/stringtie/sample_gtf/{sample}.gtf"
	log:
		"out/logs/stringtie_init/{sample}.log",
	params:
		param = config["newGene"]["stringtie_params"],
	threads: 8 
	conda:
		"../envs/stringtie.yaml"
	shell: 
		"stringtie -v -p {threads} \
		{params.param} \
		-G {input.anno} \
		-o {output} \
		{input.sbam} > {log} 2>&1"


# Merge assembled transcripts from all samples. 
# This produces a single, consistent set of transcripts and should overcome low coverage
# (leading to fragmented transcripts in some samples.
rule stringtie_merge:
	input: 
		gtf=expand("out/stringtie/sample_gtf/{sample}.gtf",sample=SAMPLES), 
		anno=config["ref"]["annotation"],
	output:
		"out/newGene/merged_assembly.gtf"
	log:
		"out/logs/stringtie_merge.log"
	threads: 8
	conda:
		"../envs/stringtie.yaml"
	shell:
		"stringtie -v --merge -p {threads} -G {input.anno} -o {output} {input.gtf} > {log} 2>&1"


# Comparison of reference annotation with all transcripts assembled by stringtie --merge.
# Thereby useful class codes are assigned describing the relation
# between reference transcripts and assembled transcripts.
rule gffcompare_transcripts:
	input:
		st_transcripts="out/newGene/merged_assembly.gtf",
		anno=config["ref"]["annotation"],
	output:
		"out/newGene/gffcompare.annotated.gtf"
	threads: 10
	conda:
		"../envs/stringtie.yaml",
	log:
		"out/logs/gffcompare_trans.log",
	script:
		"../scripts/gffcompare.py"
