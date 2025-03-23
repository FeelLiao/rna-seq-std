# Initial assembly of transcripts with StringTie for each sample 
rule stringtie_initial:
	input: 
		sbam="mapped/{sample}.bam",
		anno=config["ref"]["annotation"],
	output:
		"stringtie/sample_gtf/{sample}.gtf"
	log:
		"logs/stringtie_init/{sample}.log"
	threads: 8 
	conda:
		"../envs/stringtie.yaml"
	shell: 
		# -l adds a label to assembled transcripts. Is it necessary? 
		"stringtie -v -e -p {threads} -G {input.anno} -o {output} {input.sbam} 2> {log}"


# Merge assembled transcripts from all samples. This produces a single, consistent set of transcripts and should overcome low coverage (leading to fragmented transcripts in some samples.
rule stringtie_merge:
	input: 
		gtf=expand("stringtie/sample_gtf/{sample}.gtf",sample=SAMPLES), # or prepare a mergelist.txt of files.
		anno=config["ref"]["annotation"],
	output:
		"stringtie/stringtie_merged_assembly.gtf"
	log:
		"logs/stringtie_merge.log"
	threads: 8
	conda:
		"../envs/stringtie.yaml"
	shell:
		"stringtie -v --merge -p {threads} -G {input.anno} -o {output} {input.gtf} 2> {log}"


# Comparison of reference annotation with all transcripts assembled by stringtie --merge. Thereby usefull class codes are assigned describing the relation betwenn reference transcripts and assembled transcripts.
rule gffcompare_transcripts:
	input:
		st_transcripts="stringtie/stringtie_merged_assembly.gtf",
		anno=config["ref"]["annotation"],
	output:
		"gffcompare/GFFcompare.annotated.gtf"
	threads: 4
	conda:
		"../envs/stringtie.yaml"
	shell:
		"gffcompare -G -r {input.anno} -o gffcompare/GFFcompare {input.st_transcripts}"
