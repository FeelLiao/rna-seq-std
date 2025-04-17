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
				merged="out/quantification/samples_merged_counts.csv",
				tpm="out/quantification/samples_merged_tpm.csv",
				# fpkm="out/quantification/samples_merged_fpkm.csv"
		params:
				input_dir=lambda w, input: Path(input[0]).parent,
		log:
				"out/logs/counts_merge.log"
		conda:
				"../envs/featurecounts.yaml"
		script:
				"../scripts/merge_matrix.py"